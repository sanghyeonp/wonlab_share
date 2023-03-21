import code
# code.interact(local=dict(globals(), **locals()))

file = "c5.go.bp.v2023.1.Hs.entrez.gmt"

with open(file, 'r') as f:
    rows = [row.strip() for row in f.readlines()]

new_rows = []
for row in rows:
    ele = row.split(sep="\t")
    ele.pop(1)
    new_rows.append("\t".join(ele))

with open(file+".modified", 'w') as f:
    f.writelines([row+"\n" for row in new_rows])

# code.interact(local=dict(globals(), **locals()))
