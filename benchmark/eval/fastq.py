class LengthNotEqualException(Exception):
    pass

class Fastq():
    def __init__(self, fq=None):
        self.name = ""
        self.seq = ""
        self.plus = ""
        self.qual = ""

        if fq is not None:
            self.name = fq.name
            self.seq = fq.seq
            self.plus = fq.plus
            self.qual = fq.qual

    def get_obj(self, fp):
        self.name = fp.readline()[:-1]
        self.seq = fp.readline()[:-1]
        self.plus = fp.readline()[:-1]
        self.qual = fp.readline()[:-1]

        if self.name == "":
            return True  # eof

        #if len(self.seq) != len(self.qual):
            #raise LengthNotEqualException("Fastq sequence's length is not equal to quality score's length.")

        return False

    def __lt__(self, other):
        split1 = self.name.split("_")
        split2 = other.name.split("_")
        return int(split1[-1][:-2]) < int(split2[-1][:-2])

    def write(self, ofs):
        ofs.write(self.name + '\n')
        ofs.write(self.seq + '\n')
        ofs.write(self.plus + '\n')
        ofs.write(self.qual + '\n')