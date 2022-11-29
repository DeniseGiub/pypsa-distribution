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

#I create the bus, located in the centre of microgrid A (micA)

n.madd("Bus", ["onebus"], carrier="AC", v_nom=20)
#%%
#file "costs" last two lines are completely invended
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

def attach_conventional_generators(n, tech_modelling, conventional_carriers):

    for tech in tech_modelling:

        buses_i = n.buses.index
        n.madd(
            "Generator",
            buses_i ,
            " " + tech, 
            bus=["onebus"], 
            carrier=tech,
            p_nom_extendable=True,
            p_nom_max=100,
            marginal_cost=10, #random number (to be corrected)
            capital_cost=1000, #random number (to be corrected)
            efficiency=0.3,
            p_set=100,
            p_max_pu=1,
            )

attach_conventional_generators(
    n, 
    config["tech_modelling"]["conv_techs"],
    config["electricity"]["conventional_carriers"]
)


def attach_storageunits(n, costs,technologies, extendable_carriers ):

    for tech in technologies:
        buses_i = n.buses.index
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
load_df=pd.read_excel(r'C:\Users\denis\OneDrive\Desktop\Mini grids\pypsa-distribution\Scripts\electric_load_1.xlsx')
#%%

load_df=load_df.set_index([n.snapshots])

load_df.index.names=['time']

def attach_load(n, load_paths, load_df, tech_modelling):
    
    n.madd("Load", ["My_load"], bus=["onebus"], carrier="AC", p_set=load_df)

load_paths=r'C:\Users\denis\OneDrive\Desktop\Mini grids\pypsa-distribution\Scripts\electric_load.xlsx'

attach_load(n, load_paths, load_df, config["tech_modelling"]["load_carriers"])


print(n)

# %%
