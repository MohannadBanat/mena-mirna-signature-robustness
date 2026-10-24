import gzip, os, urllib.request
B = "https://ftp.ebi.ac.uk/pub/databases/RNAcentral/current_release"
os.makedirs("ref", exist_ok=True)
urllib.request.urlretrieve(B + "/id_mapping/database_mappings/mirbase.tsv", "ref/rc_mirbase.tsv")
urllib.request.urlretrieve(B + "/sequences/by-database/mirbase.fasta.gz", "ref/rc_mirbase.fa.gz")
m = {}
for l in open("ref/rc_mirbase.tsv"):
    f = l.rstrip("\n").split("\t")
    if len(f) >= 6 and f[3] == "9606" and f[2].startswith("MIMAT"): m.setdefault(f[0], []).append((f[5], f[2]))
seq = {}; cur = None
for l in gzip.open("ref/rc_mirbase.fa.gz", "rt"):
    if l.startswith(">"):
        cur = l[1:].split()[0].split("_")[0]
        if cur in m: seq[cur] = []
    elif cur in seq: seq[cur].append(l.strip())
n = 0
with open("ref/mature.fa", "w") as o:
    for u, lst in m.items():
        if not seq.get(u): continue
        s = "".join(seq[u]).replace("T", "U")
        for name, acc in sorted(set(lst)):
            o.write(">%s %s Homo sapiens %s\n%s\n" % (name, acc, name[4:], s)); n += 1
print("RNAcentral hsa mature written:", n)
