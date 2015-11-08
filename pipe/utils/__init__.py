import csv
from clint.textui import puts_err, colored, indent
import os


def check_file(filename):
    # Checks whether a file exists
    if not os.path.exists(filename):
        with indent(4):
            exit(puts_err(colored.red("\n" + filename + "not found\n")))
    return filename

def ensure_dir(directory):
    if not os.path.exists(directory):
        os.makedirs(directory)
    return directory

class sample_file:
    """
        The sample file holds all the information
        regarding every sequence run and 
        associated groupings.
    """
    def __init__(self, sample_file, config):
        self.filename = check_file(sample_file)
        self.file_in = open(sample_file, 'rU')
        self.records = []
        self.config = config

        if config["debug"]:
            self.config["analysis_path"] = config["analysis_path"] +"_debug"

        # Prepend directories
        self.config["bam_path"] = ensure_dir(self.config["analysis_path"] + "/" + self.config["bam_path"])
        self.config["vcf_path"] = ensure_dir(self.config["analysis_path"] + "/" + self.config["vcf_path"])
        self.config["results_path"] = ensure_dir(self.config["analysis_path"] + "/" + self.config["results_path"])

        # Generate any additional necessary directories
        ensure_dir(self.config["results_path"] + "/dups")
        ensure_dir(self.config["results_path"] + "/telseq")
        ensure_dir(self.config["results_path"] + "/jellyfish")
        ensure_dir(self.config["results_path"] + "/eav")
        ensure_dir(self.config["vcf_path"] + "/sv")
        ensure_dir(self.config["vcf_path"] + "/snps")
        ensure_dir(self.config["tmp_path"] + "/dups")

        dialect = csv.Sniffer().sniff(self.file_in.read(10000))
        self.file_in.seek(0)
        for record in csv.DictReader(self.file_in, dialect = dialect):
            # Setup read group
            rg = [r"{k}:{v}".format(k=x, v=record[x]) for x in ["ID","SM","PL"] if record[x] != ""]

            record["readgroup"] = r"@RG\t" + r'\t'.join(rg)
            record["FQ1"] = self.config["fastq_path"] + "/" + record["FQ1"]
            record["FQ2"] = self.config["fastq_path"] + "/" + record["FQ2"]
            if config["debug"]:
                record["FQ1"] = "<( zcat " + record["FQ1"] + " | head -n " + str(config["debug_reads"]*4) + " )"
                record["FQ2"] = "<( zcat " + record["FQ2"] + " | head -n " + str(config["debug_reads"]*4) + " )"
            record["bam_individual"] = self.config["bam_path"] + "/" + record["ID"] + ".bam"
            self.records.append(record)

            # Check that IDs are unique
            IDs = [x["ID"] for x in self.records]
            try:
                assert(len(IDs) == len(set(IDs)))
            except:
                with indent(4):
                    exit(puts_err(colored.red("\nIDs in sample set are not unique\n")))

        if config["debug"]:
            self.records = self.records[0:config["debug_samples"]]
        

        self.sample_groups = list(set([x["SM"] for x in self.records]))

    def iterate_fq(self):
        """
            Iterate through fastq pairs
        """
        for i in self.records:
            yield {"FQ1": self.config["fastq_path"] + "/" + i["FQ1"], 
                   "FQ2": self.config["fastq_path"] + "/" + i["FQ2"],
                   "readgroup" : i["readgroup"],
                   "bam_individual": self.config["bam_path"] + "/" + i["readgroup"]}

    def iterate_SM_sets(self):
        for SM in self.sample_groups:
            SM_set = {}
            SM_set["SM"] = SM
            SM_set["fq_set"] = [x for x in self.records if x["SM"] == SM]
            SM_set["bamlist"] = ' '.join([x["bam_individual"] for x in self.records if x["SM"] == SM])
            SM_set["bam_merged"] = self.config["bam_path"] + "/" + SM + ".bam"
            yield SM_set

    def iterate_bams(self):
        for SM in self.sample_groups:
            yield self.config["bam_path"] + "/" + SM + ".bam"

    def iterate_vcfs(self):
        for SM in self.sample_groups:
            yield self.config["vcf_path"] + "/" + SM + ".vcf.gz"

