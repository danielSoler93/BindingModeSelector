import os
import glob
import filecmp
from BMSelector import main as mn


TEST_DIR = os.path.dirname(os.path.abspath(__file__))
PDB = glob.glob(os.path.join(TEST_DIR, "data/pdb/*.pdb"))[0]
CPPTRAJ = glob.glob(os.path.join(TEST_DIR, "data/cpptraj/*"))[0]
TLEAP = glob.glob(os.path.join(TEST_DIR, "data/tleap/*"))[0]
FP = os.path.join(TEST_DIR, "data/fp/")
COMPARE = os.path.join(TEST_DIR, "data/compare")


def check_cpptraj(file1, file2, folder):
    with open(file1, 'r') as f1:
        lines1 = f1.readlines()
    with open(file2, 'r') as f2:
        lines2 = f2.readlines()
    for l1, l2 in zip(lines1, lines2):
        if l1 == l2:
           assert True
        else:
           residue_part1, output_part1 = l1.split("out ")
           residue_part2, output_part2 = l2.split("out ")

           _, r11, r12 = residue_part1.split(":")
           _, r21, r22 = residue_part2.split(":")
           assert int(r11.strip()) == int(r21.strip())-1
           assert int(r12.strip()) == int(r22.strip())-1

           assert output_part2.strip() == os.path.join(folder, output_part1.strip())
           
def check_tleap(file1, file2, pdb_in, pdb_out, folder):
    with open(file1, 'r') as f1:
        lines1 = f1.readlines()
    with open(file2, 'r') as f2:
        lines2 = f2.readlines()
    for l1, l2 in zip(lines1, lines2):
        if l1 == l2:
           assert True
        else:
           if "loadpdb" in l1:
               pdb1 = l1.split()[-1]
               pdb2 = l2.split()[-1]
               assert pdb1.strip() == "$PDB"
               assert pdb2.strip() == pdb_in
           elif "complex" in  l1:
               pdb1 = l1.split()[-1]
               pdb2 = l2.split()[-1]
               assert pdb1.strip() == "$OUTPUT"
               assert pdb2.strip() == pdb_out
           

def test_pdb_clean_name(pdb=PDB):
    out_file, folder = mn.pdb_clean_name(pdb)
    condition1 = "-" not in out_file 
    condition2 = "." not in out_file[:-4] 
    condition3 = os.path.exists(out_file)
    condition4 = filecmp.cmp(pdb, out_file)
    condition5 = folder == os.path.basename(out_file)[:-4]
    assert condition1 
    assert condition2 
    assert condition3 
    assert condition4
    assert condition5


def test_reoder_cpptraj_input(cpptraj=CPPTRAJ, folder="output"):
    out_file = mn.reoder_cpptraj_input(cpptraj, folder)
    check_cpptraj(cpptraj, out_file, folder)
    
def test_template_tleap(template=TLEAP, pdb_in="a.pdb", pdb_out="b.pdb", folder_out="output"):
    tleap_input = mn.template_tleap(template, pdb_in, pdb_out, folder_out)
    check_tleap(template, tleap_input, pdb_in, pdb_out, folder_out)

def test_get_fp(folder=FP):
    _, vdw, elec = mn.get_fp(folder)[0]
    assert vdw == -0.2541
    assert elec == -12.4371

def test_compare_fp(folder=FP, compare_folder=COMPARE):
    fp = mn.get_fp(folder)
    comparation, extrems = mn.compare_fp(fp, compare_folder)
    assert extrems == [-0.6046547948110378, 0.6034926588644364, -12.7548549863448, -6.955548838559576]
    assert comparation[0] == True
    
def test_main(pdbs=PDB):
    results = mn.main(pdbs, CPPTRAJ, TLEAP, COMPARE)
    assert results[PDB] == 5
