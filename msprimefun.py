import stdpopsim
species = stdpopsim.get_species("HomSap")
#import the species map homosapien
#will now attempt to modify "isolationwithmigration to include 3 populations instead of 2.
#still need to work out how to do the split
import copy
import textwrap
import msprime
class DemographicModel:
   
    def __init__(
        self,
        *,
        id,
        description,
        long_description,
        generation_time=None,
        citations=None,
        qc_model=None,
        population_configurations=None,
        demographic_events=None,
        migration_matrix=None,
        population_id_map=None,
        populations=None,
        model=None,
        mutation_rate=None,
    ):
        self.id = id
        self.description = description
        self.long_description = long_description
        self.generation_time = 1 if generation_time is None else generation_time
        self.mutation_rate = mutation_rate
        self.citations = [] if citations is None else citations
        self.qc_model = qc_model

        if model is None:
            assert population_configurations is not None
            assert populations is not None
            assert len(populations) == len(population_configurations)
            population_configurations = copy.deepcopy(population_configurations)
            
            for pop, pop_config in zip(populations, population_configurations):
                if pop_config.metadata is None:
                    pop_config.metadata = {}
                if pop.id != "" and not pop.id.startswith("qc_"):
                    pop_config.metadata["name"] = pop.id
                    pop_config.metadata["description"] = pop.description
            
            self.model = msprime.Demography.from_old_style(
                population_configurations=population_configurations,
                demographic_events=demographic_events,
                migration_matrix=migration_matrix,
                population_map=population_id_map,
            )
            for msp_pop, local_pop in zip(self.model.populations, populations):
                
                msp_pop.allow_samples = True
                if local_pop.sampling_time is not None:
                    msp_pop.default_sampling_time = local_pop.sampling_time
                else:
                    msp_pop.allow_samples = False
        else:
            assert population_configurations is None
            assert population_id_map is None
            assert populations is None
            self.model = copy.deepcopy(model)
            # See note above. We allow samples in all populations for now.
            for pop in self.model.populations:
                pop.allow_samples = True

    def __str__(self):
        long_desc_lines = [
            line.strip()
            for line in textwrap.wrap(textwrap.dedent(self.long_description))
        ]
        long_desc = "\n║                     ".join(long_desc_lines)
        s = (
            "Demographic model:\n"
            f"║  id               = {self.id}\n"
            f"║  description      = {self.description}\n"
            f"║  long_description = {long_desc}\n"
            f"║  generation_time  = {self.generation_time}\n"
            f"║  mutation_rate    = {self.mutation_rate}\n"
            f"║  citations        = {[cite.doi for cite in self.citations]}\n"
            f"║{self.model}"
        )
        return s

    @property
    def populations(self):
        return self.model.populations

    @property
    def num_populations(self):
        return self.model.num_populations

    @property
    def num_sampling_populations(self):
        return sum(int(pop.allow_samples) for pop in self.populations)


    def register_qc(self, qc_model):
        if not isinstance(qc_model, self.__class__):
            raise ValueError(f"Cannot register non-DemographicModel '{qc_model}' as QC mode.")
        if self.qc_model is not None:
            raise ValueError(f"QC model already registered for {self.id}.")
        self.qc_model = qc_model

    def get_samples(self, *args):
        samples = []
        for pop_index, n in enumerate(args):
            if self.populations[pop_index].allow_samples:
                samples.append(
                    msprime.SampleSet(
                        num_samples=n,
                        population=pop_index,
                        time=self.populations[pop_index].default_sampling_time,
                        ploidy=1,  # Avoid breaking too much at once.
                    )
                )
            elif n > 0:
                raise ValueError(
                    "Samples requested from non-sampling population {pop_index}"
                )
        return samples
#name __init__ is not defined??
class isolationwithmigration:
    def __init__(self, NA, N1, N2, N3, T, M12, M21, M13, M31, M23, M32):
        model=msprime.Demography()
        model.add_population(initial_size=N1, name="pop1")
        model.add_population(initial_size=N2, name="pop2")
        model.add_population(initial_size=N3, name="pop3")
        model.add_population(initial_size=NA, name="ancestral")
        model.set_migration_rate(source="pop1", dest="pop2", rate=M12)
        model.set_migration_rate(source="pop2", dest="pop1", rate=M21)
        model.add_population_split(
            time=T, ancestral="ancestral", derived=["pop1", "pop2","pop3"]
        )
        long_description = ""
#maybe add more splits here?
        super().__init__(
        id="isolationwithmigration",
        description="Generic IM model",
        long_description=long_description,
        model=model,
        generation_time=1,
    )
#end of trying to configure new demographic model
#below is 3 populations, need to figure out what to use for second split
model = stdpopsim.isolationwithmigration(NA=7500, N1=2500, N2=2500, N3=2500, T=10, M12=0, M21=0, M13=0, M31=0, M23=0, M32=0)
contig = species.get_contig("chr22")
samples=model.get_samples(10,10,10)
engine=stdpopsim.get_engine("msprime")
ts=engine.simulate(model, contig, samples, dryrun=True)

#3 new errors!
#when run, says "stdpopsim" has not attribute 'isolationwithmigration
#when replaceing "stdpopsim" with "species", says "species" has no attribute "isolationwithmigration"
#if deleting species, says "TypeError: object.__init__() takes exactly one argument"
