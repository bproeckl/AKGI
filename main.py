aminoDict = {"TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
             "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
             "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
             "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
             "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
             "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
             "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
             "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
             "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
             "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
             "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
             "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
             "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
             "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
             "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
             "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G"}
endcodons = ["TAA", "TAG", "TGA"]


class Read:
    def __init__(self, lines):
        self.name = lines[0][1:]
        self.bases = lines[1]
        self.bases_complementary = self.bases[::-1]
        self.bases_complementary = self.bases_complementary.replace("T", "U").replace("A", "T").replace("U", "A")
        self.bases_complementary = self.bases_complementary.replace("G", "U").replace("C", "G").replace("U", "C")

    def get_atgs(self):
        atgs = []
        atgscomplementary = []
        for x in range(0, len(self.bases) - 3):
            if self.bases[x:x + 3] == "ATG":
                atgs.append(x)
            if self.bases_complementary[x:x + 3] == "ATG":
                atgscomplementary.append(x)
        return atgs, atgscomplementary

    def __str__(self):
        return self.name + ": " + self.bases + ": " + self.bases_complementary

    def __repr__(self):
        return self.__str__()

    def __eq__(self, other):
        if not isinstance(other, Read):
            return False
        if other.name == self.name and other.bases == self.bases:
            return True
        return False


class ORF:
    def __init__(self, codons):
        self.dna=""
        self.aa=""
        for x in codons:
            self.dna+=x
            self.aa+=aminoDict[x]

    def __str__(self):
        return self.dna + ": " + self.aa

    def __repr__(self):
        return self.__str__()

class Decoder:
    def __init__(self, read, minlen, maxlen):
        self.read = read
        self.orfs = {}
        self.minlen = minlen
        self.maxlen = maxlen

    def get_orfs(self):
        atgs, atgscomplementary = self.read.get_atgs()
        for x in atgs:
            codons=self.get_codons(x)
            if codons!=[]:
                self.orfs[x+1]=ORF(codons)
        for x in atgscomplementary:
            codons = self.get_codons(x,True)
            if codons != []:
                self.orfs[-len(self.read.bases)+x] = ORF(codons)
        print(self.orfs)


    def get_codons(self, start, complementary=False):
        if self.minlen * 3 + start + 3 > len(self.read.bases):
            return []
        codons = []
        for i in range(0, self.maxlen - 1):
            x = i * 3 + start
            if x < len(self.read.bases):
                if complementary:
                    codon = self.read.bases_complementary[x:x + 3]
                else:
                    codon = self.read.bases[x:x + 3]
                codons.append(codon)
                if codon in endcodons:
                    break
        if len(codons)<=self.minlen:
            return []
        else:
            return codons
    def write_dna(outfile_dna):
        dna=""
        f = open(outfile_dna, "w")
        f.write(dna)
        f.close()

    def write_aa(outfile_aa):
        aa = ""
        f = open(outfile_aa, "w")
        f.write(aa)
        f.close()


def read_fasta(readfile):
    file = open(readfile)
    data = file.read()
    data = data.split(">")
    data.pop(0)
    reads = []
    for x in data:
        x = x.split("\n")
        x[0] = ">" + x[0]
        reads.append(Read(x))
    return reads


# def findORFs(infile, outfile_dna, outfile_aa, minlen, maxlen):
    # writefasta

if __name__ == "__main__":
    reads = read_fasta("test.fasta")
    print(reads[0])
    print(reads[0].get_atgs())
    dec = Decoder(reads[0], 1, 10)
    dec.get_orfs()
