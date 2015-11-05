import csv
from clint.textui import puts_err, colored, indent

class sample_file:
    """
        The sample file holds all the information
        regarding every sequence run and 
        associated groupings.
    """
    def __init__(self, sample_file, config):
        self.filename = sample_file
        self.file_in = open(sample_file, 'rU')
        self.records = []
        self.config = config
        dialect = csv.Sniffer().sniff(self.file_in.read(1024))
        self.file_in.seek(0)
        for record in csv.DictReader(self.file_in, dialect = dialect):
            rg = [r"{k}:{v}".format(k=x, v=record[x]) for x in ["ID","SM","PL"] if record[x] != ""]
            record["readgroup"] = r"@RG\t" + r'\t'.join(rg)
            self.records.append(record)

            # Check that IDs are unique
            IDs = [x["ID"] for x in self.records]
            try:
                assert(len(IDs) == len(set(IDs)))
            except:
                with indent(4):
                    exit(puts_err(colored.red("\nIDs in sample set are not unique\n")))

        print self.records
        self.sample_groups = list(set([x["SM"] for x in self.records]))

    def iterate_fq(self):
        """
            Iterate through fastq pairs
        """
        for i in self.records:
            yield {"fq1": self.config["fastq_path"] + "/" + i["FQ1"], 
                   "fq2": self.config["fastq_path"] + "/" + i["FQ2"]}

    def iterate_SM_sets(self):
        for SM in self.sample_groups:
            yield {SM: [x for x in self.records if x["SM"] == SM]}

    def iterate_bams(self):
        for SM in self.sample_groups:
            yield self.config["bam_path"] + "/" + SM + ".bam"

    def iterate_vcfs(self):
        for SM in self.sample_groups:
            yield self.config["vcf_path"] + "/" + SM + ".vcf.gz"

