
import subprocess
from pymol import cmd as pm

resi_number = {
    '2TPI':(1016,'Z'),
    
}

for obj, (resi, chain) in resi_number.items():
    pm.fetch(obj)
    pm.alter(f'bymolecule {obj} and resi {resi} and chain {chain}', 'chain="P"')
    arquivo_pdb = f'/home/gessualdo/gessualdo/pymol_pdb/{obj}.pdb'
    pm.save(arquivo_pdb, obj)

    print(obj, chain)
    caminho = '/home/gessualdo/anaconda3/envs/envname/bin/plip'
    plip = subprocess.run([caminho, '-f', arquivo_pdb, '--peptide', 'P',
                           '-o', f'/home/gessualdo/gessualdo/plip_pdbs/{obj}', '-x'])



