import msprime
demography = msprime.Demography()
demography.add_population(name="NA1", initial_size=10_000)
demography.add_population(name="NA2", initial_size=5_000)
demography.add_population(name="N1", initial_size=1_000)
demography.add_population(name="N2", initial_size=10_000)
demography.add_population(name="N3", initial_size=5_000)
demography.add_population_split(time=1000, derived=["NA2", "N1"], ancestral="NA1")
demography.add_population_split(time=500, derived=["N2", "N3"], ancestral="NA2")
demography.set_migration_rate(source="N1", dest="N2", rate=0.1)
demography.set_migration_rate(source="N2", dest="N1", rate=0.1)
demography.set_migration_rate(source="N1", dest="N3", rate=0.1)
demography.set_migration_rate(source="N3", dest="N1", rate=0.1)
demography.set_migration_rate(source="N2", dest="N3", rate=0.1)
demography.set_migration_rate(source="N3", dest="N2", rate=0.1)
demography.sort_events()
ts = msprime.sim_ancestry(samples={"N1": 1, "N2": 1, "N3":1}, demography=demography, random_seed=12)
ts

print(demography)
import demesdraw
import matplotlib.pyplot as plt 
graph = msprime.Demography.to_demes(demography)
fig, ax = plt.subplots()  # use plt.rcParams["figure.figsize"]
demesdraw.tubes(graph, ax=ax, seed=1)
plt.show() 



