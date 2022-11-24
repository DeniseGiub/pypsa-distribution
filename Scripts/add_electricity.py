#Adding components to the network 
#%%
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import yaml
from create_network import n

#Load yaml file
with open(f"C://Users//denis//OneDrive//Desktop//Mini grids//pypsa-distribution//Scripts//config.yaml") as file:
        config = yaml.load(file, Loader=yaml.FullLoader)
#%%
#I create the bus, located in the centre of microgrid A (micA)

n.madd("Bus", ["onebus"], x=10.2, y=9.3, carrier="AC", v_nom=20)
# x=config["microgrids_list"]["micA"]["lon"], y=config["microgrids_list"]["micA"]["lat"]

#TODO file cost last two lines are completely invended
costs=pd.read_csv(r'C:\Users\denis\OneDrive\Desktop\Mini grids\pypsa-distribution\costs.csv')

tech_costs=pd.read_csv(r'C:\Users\denis\OneDrive\Desktop\Mini grids\pypsa-distribution\costs.csv')

idx = pd.IndexSlice

Nyears = n.snapshot_weightings.objective.sum() / 8760.0

def calculate_annuity(n, r):
    """
    Calculate the annuity factor for an asset with lifetime n years and
    discount rate of r, e.g. annuity(20, 0.05) * 20 = 1.6
    """
    if isinstance(r, pd.Series):
        return pd.Series(1 / n, index=r.index).where(
            r == 0, r / (1.0 - 1.0 / (1.0 + r) ** n)
        )
    elif r > 0:
        return r / (1.0 - 1.0 / (1.0 + r) ** n)
    else:
        return 1 / n

def _add_missing_carriers_from_costs(n, costs, carriers):
    missing_carriers = pd.Index(carriers).difference(n.carriers.index)
    if missing_carriers.empty:
        return

    emissions_cols = (
        costs.columns.to_series().loc[lambda s: s.str.endswith("_emissions")].values
    )
    suptechs = missing_carriers.str.split("-").str[0]
    emissions = costs.loc[suptechs, emissions_cols].fillna(0.0)
    emissions.index = missing_carriers

tech_costs='C://Users//denis//OneDrive//Desktop//Mini grids//pypsa-distribution//costs.csv'

def load_costs(tech_costs, config, elec_config, Nyears=1):
    """
    set all asset costs and other parameters
    """
    costs = pd.read_csv(tech_costs, index_col=list(range(3))).sort_index()

    # correct units to MW and EUR
    costs.loc[costs.unit.str.contains("/kW"), "value"] *= 1e3
    costs.loc[costs.unit.str.contains("USD"), "value"] *= config["USD2013_to_EUR2013"]

    costs = (
        costs.loc[idx[:, config["year"], :], "value"]
        .unstack(level=2)
        .groupby("technology")
        .sum(min_count=1)
    )

    costs = costs.fillna(
        {
            "CO2 intensity": 0,
            "FOM": 0,
            "VOM": 0,
            "discount rate": config["discountrate"],
            "efficiency": 1,
            "fuel": 0,
            "investment": 0,
            "lifetime": 25,
        }
    )

    costs["capital_cost"] = (
        (
            calculate_annuity(costs["lifetime"], costs["discount rate"])
            + costs["FOM"] / 100.0
        )
        * costs["investment"]
        * Nyears
    )

    costs.at["OCGT", "fuel"] = costs.at["gas", "fuel"]
    costs.at["CCGT", "fuel"] = costs.at["gas", "fuel"]

    costs["marginal_cost"] = costs["VOM"] + costs["fuel"] / costs["efficiency"]

    costs = costs.rename(columns={"CO2 intensity": "co2_emissions"})

    costs.at["OCGT", "co2_emissions"] = costs.at["gas", "co2_emissions"]
    costs.at["CCGT", "co2_emissions"] = costs.at["gas", "co2_emissions"]

    costs.at["solar", "capital_cost"] = 0.5 * (
        costs.at["solar-rooftop", "capital_cost"]
        + costs.at["solar-utility", "capital_cost"]
    )

    def costs_for_storage(store, link1, link2=None, max_hours=1.0):
        capital_cost = link1["capital_cost"] + max_hours * store["capital_cost"]
        if link2 is not None:
            capital_cost += link2["capital_cost"]
        return pd.Series(
            dict(capital_cost=capital_cost, marginal_cost=0.0, co2_emissions=0.0)
        )
    
    max_hours = elec_config["max_hours"]           
    costs.loc["battery"] = costs_for_storage(                      
            costs.loc["lithium"],                            #line 119 in file costs.csv' which was battery storage was modified into lithium (same values left)
            costs.loc["battery inverter"],
            max_hours=max_hours["battery"],
    )
    max_hours = elec_config["max_hours"]
    costs.loc["battery"] = costs_for_storage(
            costs.loc["lead acid"],                          #line 120 in file 'costs.csv' which was battery storage was modified into lithium (same values left)
            costs.loc["battery inverter"],
            max_hours=max_hours["battery"],
    )
    
    costs.loc["H2"] = costs_for_storage(
        costs.loc["hydrogen storage"],
        costs.loc["fuel cell"],
        costs.loc["electrolysis"],
        max_hours=max_hours["H2"],
    )

    for attr in ("marginal_cost", "capital_cost"):
        overwrites = config.get(attr)
        if overwrites is not None:
            overwrites = pd.Series(overwrites)
            costs.loc[overwrites.index, attr] = overwrites

    return costs

costs = load_costs(
    tech_costs,
    config["costs"],
    config["electricity"],
    Nyears,
    )
#%%
def attach_wind_and_solar(n, costs, tech_modelling, extendable_carriers):

    _add_missing_carriers_from_costs(n, costs, tech_modelling)

    for tech in tech_modelling:
       
        with xr.open_dataset(f"renewable_profiles/profile_{tech}.nc") as ds:
            
            if ds.indexes["bus"].empty:
                continue   

            suptech = tech.split("-", 2)[0]
          
            n.madd(
            "Generator",
            ds.indexes["bus"],
            " " + tech,
            bus=["onebus"],
            carrier=tech,
            p_nom_extendable=tech in extendable_carriers["Generator"],
            p_nom_max=ds["p_nom_max"].to_pandas(), #look at the config 
            weight=ds["weight"].to_pandas(),
            marginal_cost=costs.at[suptech, "marginal_cost"],
            capital_cost=costs.at[tech, "capital_cost"],
            efficiency=costs.at[suptech, "efficiency"],
            p_set=ds["profile"].transpose("time", "bus").to_pandas().reindex(n.snapshots),
            p_max_pu=ds["profile"].transpose("time", "bus").to_pandas().reindex(n.snapshots),
            )

attach_wind_and_solar(
    n,
    costs,
    config["tech_modelling"]["general_vre"],
    config["electricity"]["extendable_carriers"],
    )
#%%
def load_powerplants(ppl_fn):
    carrier_dict = {
        "ocgt": "OCGT",
        "ccgt": "CCGT",
        "bioenergy": "biomass",
        "ccgt, thermal": "CCGT",
        "hard coal": "coal",
    }
    return (
        pd.read_csv(ppl_fn, index_col=0, dtype={"bus": "str"})
        # .to_pypsa_names() #commented because I had error "AttributeError: 'DataFrame' object has no attribute 'to_pypsa_names'"
        # .convert_country_to_alpha2() #commented because I had error "AttributeError: 'DataFrame' object has no attribute 'country_to_Alpha2'"
        .rename(columns=str.lower)
        .drop(columns=["efficiency"])
        .replace({"carrier": carrier_dict})
    )


ppl_fn="C://Users//denis//OneDrive//Desktop//Mini grids//pypsa-distribution//Scripts//powerplants.csv"

ppl=load_powerplants(ppl_fn)
#%%

# def attach_conventional_generators(
#     n,
#     costs,
#     ppl,
#     tech_modelling,
#     extendable_carriers,
#     conventional_carriers,
# ):
#     carriers={'oil'}
#     _add_missing_carriers_from_costs(n, costs, tech_modelling)

#     ppl = (
#         ppl.query("fueltype in @carriers")
#         .join(costs, on="fueltype", rsuffix="_r")
#         .rename(index=lambda s: "C" + str(s))
#     )
#     ppl["efficiency"] = ppl.efficiency.fillna(ppl.efficiency)
    
#     n.madd(
#         "Generator",
#         ppl.index,
#         carrier=ppl.technology,
#         bus=["onebus"],
        # p_nom_min=ppl.p_nom.where(ppl.rechnology.isin(conventional_carriers), 0),
        # p_nom=ppl.p_nom.where(ppl.technology.isin(conventional_carriers), 0),
        # p_nom_extendable=ppl.technology.isin(extendable_carriers["Generator"]),
        # efficiency=ppl.efficiency,
        # marginal_cost=ppl.marginal_cost,
        # capital_cost=ppl.capital_cost,
        # build_year=ppl.datein.fillna(0).astype(int),
        # lifetime=(ppl.dateout - ppl.datein).fillna(np.inf),

    
   
    
#%%
# attach_conventional_generators(
#     n,
#     costs,
#     ppl,
#     config["tech_modelling"]["conv_techs"],
#     config["electricity"]["extendable_carriers"],
#     config["electricity"]["conventional_carriers"],
#     )
    
#%%
def attach_storageunits(n, costs,technologies, extendable_carriers ):

    buses_i = n.buses.index

    for tech in technologies:
        
        n.madd(
            "StorageUnit",
            buses_i, 
            " " + tech, 
            bus=["onebus"],
            carrier=tech,
            p_nom_extendable=True,
            capital_cost=costs.at[tech, "capital_cost"],
            marginal_cost=costs.at[tech, "marginal_cost"],
            # efficiency_store=costs.at[lookup_store[tech], "efficiency"],
            # efficiency_dispatch=costs.at[lookup_dispatch[tech], "efficiency"],
            # max_hours=max_hours[tech],
            cyclic_state_of_charge=True
            
        )

attach_storageunits(n, 
                    costs,
                    config["tech_modelling"]["storage_techs"],
                    config["electricity"]["extendable_carriers"],
                    )
#%%
load_df=pd.read_excel(r'C:\Users\denis\OneDrive\Desktop\Mini grids\pypsa-distribution\Scripts\electric_load.xlsx')
#%%
load_df=load_df(header=None)
#%%
# n_load=35
load_df=load_df.set_index([n.snapshots])
#%%
load_df.index.names=['time']
#%%
# index=pd.Index( list(range(n_load)))

def attach_load(n, load_paths, load_df, tech_modelling):
    
    n.madd("Load", ["My_load"], bus=["onebus"], carrier="AC", p_set=load_df)

load_paths=r'C:\Users\denis\OneDrive\Desktop\Mini grids\pypsa-distribution\Scripts\electric_load.xlsx'

attach_load(n, load_paths, load_df, config["tech_modelling"]["load_carriers"])

print(n)
#%%
# Optimization
from pypsa.linopf import ilopf

solver_name="gurobi"

n.lopf(n.snapshots, solver_name=solver_name, pyomo=False)

# %%
