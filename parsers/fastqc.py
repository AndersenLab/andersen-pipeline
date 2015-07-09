from db.eav import eav

class fastqc(object):
    def __init__(self, filename, sample = "", analysis = ""):
        self.filename = filename
        self.sample = sample
        self.analysis = analysis
        self.summary = {}
        self.basic_statistics = {}
        with open(self.filename) as f:
            for line in f:
                line = line.strip()
                if line.startswith("##"):
                    pass
                elif line.find("END_MODULE") > 0:
                    pass
                elif line.startswith(">>"):
                    module,v = line.strip(">").split("\t")
                    self.summary[module] = v
                elif line.startswith("#"):
                    cols = line.strip("#").split("\t")
                else:
                    if module == "Basic Statistics":
                        key, val = line.split("\t")
                        self.basic_statistics[key] = val
                print(line)

    def to_eav(self):
        # Parse summary
        for k,v in self.summary.items():
            eav(entity = self.sample, 
                attribute = k,
                value = v,
                program = "fastqc").save()
        for k,v in self.basic_statistics.items():
            print(k,v)
            eav(entity = self.sample, 
                attribute = k,
                value = v,
                program = "fastqc").save()


#x = fastqc("../results/fastqc/BGI2-RET2-test1-paired-1.fq_fastqc/fastqc_data.txt", "Test")

#print x.summary

#print(x.to_eav())