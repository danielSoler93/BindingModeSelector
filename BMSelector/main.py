import os
from multiprocessing import Pool
import time
from tqdm import tqdm 
import numpy as np
import pandas as pd
import shutil
import glob
import string



def reorder(inputfile, outputfile):
    new_pdb = os.path.basename(inputfile)[:-4] + "_reorder.pdb"
    new_pdb = new_pdb.replace("-", "_")
    folder_name = new_pdb[:-4]
    if not os.path.exists(folder_name): os.mkdir(folder_name)
    new_pdb = os.path.join(folder_name, new_pdb)

    with open(inputfile, "r") as f:
        residues = list(set([int(line[22:27].strip()) for line in f if (line.startswith("ATO") or line.startswith("HET"))]))
    
    with open(outputfile, "r") as f:
        new_lines = []
    
        current_res = "FDRDA"
        res = -1
        for line in f:
            if line.startswith("ATO") or line.startswith("HET"):
                if line[22:27].strip() != current_res:
                    current_res = line[22:27].strip()
                    res += 1  
                new_line = line[:23] + str(residues[res]).rjust(3, " ") + line[26:]
                
            else:
                new_line = line
            new_lines.append(new_line)
    
    
    with open(new_pdb, "w") as f:
        f.write(("".join(new_lines)))
    return new_pdb, folder_name


def get_fp(folder_name):
    final_fps = []
    fps = glob.glob(os.path.join(folder_name, "fp_*.txt"))
    for fp in fps:
        df=pd.read_table(fp, sep='\s+',names=['Frames','VDW', 'ELEC'], skiprows=1)
        means_vdw = np.mean(df["VDW"].values)
        means_elec = np.mean(df["ELEC"].values)
        info = [fp, means_vdw, means_elec]
        final_fps.append(info) 
    return final_fps

def template_cpptraj(cpptraj_input, new_pdb, folder):
    with open(cpptraj_input, "r") as f:
        content = string.Template("".join(f.readlines()))
        content = content.substitute(PDB=new_pdb)
    out_file = os.path.join(folder, cpptraj_input)
    with open(out_file, "w") as f:
        f.write(content)
    return out_file

def template_tleap(tleap_input, new_pdb, output, folder):
    if not os.path.exists(folder): os.mkdir(folder)
    with open(tleap_input, "r") as f:
        content = string.Template("".join(f.readlines()))
        content = content.substitute(PDB=new_pdb, OUTPUT=output)
    out_file = os.path.join(folder, os.path.basename(tleap_input))
    with open(out_file, "w") as f:
        f.write(content)
    return out_file


def compare_fp(fp, fp_compare_folder="fp_BU_MD/"):
    comparation = []
    folder = fp_compare_folder
    files = np.load(os.path.join(folder, "files.npy"))
    vdw_means = np.load(os.path.join(folder, "means_vdw.npy"))
    vdw_stds = np.load(os.path.join(folder, "std_vdw.npy"))
    elec_means = np.load(os.path.join(folder, "means_elec.npy"))
    elec_stds = np.load(os.path.join(folder, "std_elec.npy"))
    for filename, vdw, elec in fp:
        for filename_BU, vdw_mean, vdw_std, elec_mean, elec_std in zip(files, vdw_means, vdw_stds, elec_means, elec_stds):
            if os.path.basename(filename) == os.path.basename(filename_BU):
                extreme_vdw1 = vdw_mean - vdw_std
                extreme_vdw2 = vdw_mean + vdw_std
                extreme_elec1 = elec_mean - elec_std
                extreme_elec2 = elec_mean + elec_std
                if vdw > extreme_vdw1 and vdw < extreme_vdw2 and elec > extreme_elec1 and elec<extreme_elec2:
                    comparation.append(True)
                else:
                    comparation.append(False)
    return comparation, [extreme_vdw1, extreme_vdw2, extreme_elec1, extreme_elec2]
        
      
def reoder_cpptraj_input(cpptrajinput, folder, gap=1):
    """
    There is a 1 residue gap between MD and PELE 
    apply that to cpptraj input
    """
    output = os.path.basename(cpptrajinput).split(".")[0] + "_reorder.txt"
    newlines = []
    with open(cpptrajinput, "r") as f:
        lines = f.readlines()
    for l in lines:
         info = l.split("out")
         if len(info) > 1:
             first_sentence, second_sentence = info
             part1, part2, part3 = first_sentence.split(":")
             old_res = [part2, part3]
             new_res = [str(int(r.strip()) + gap) for r in old_res]
             new_residues = [o.replace(o.strip(), n) for o, n in zip(old_res, new_res)]
             first_sentence = ":".join([part1, new_residues[0], new_residues[1]])
             newline = "out".join([first_sentence, " " + os.path.join(folder, second_sentence.strip()) + "\n"])
             newlines.append(newline)
         else:
             newlines.append(l)
    with open(output, "w") as f:
        f.write("".join(newlines))
    return output
        
        
def fix_pdb(folder, pdb):
    outpdb = os.path.join(folder, os.path.basename(pdb)[:-4] + "clean.pdb")
    newlines = []
    with open(pdb, "r") as f:
        for l in f:
            if "ATOM      9  H   MET R  65" in l:
                pass
            else:
                newlines.append(l)
    with open(outpdb, "w") as f:
        f.write("".join(newlines))
    return outpdb


def pdb_clean_name(pdb):
    filename = os.path.basename(pdb)
    filename_clean = filename[:-4].replace("-", "_").replace(".", "")+".pdb"
    folder = filename_clean[:-4]
    if not os.path.exists(folder): os.mkdir(folder)
    output_file = os.path.join(folder, filename_clean)
    shutil.copy(pdb, output_file)
    return output_file, folder

def extract_info_file(filename, cpptraj_template, leap_template, MD_folder, gap):
        filename, folder = pdb_clean_name(filename)
        cpptraj_clean = reoder_cpptraj_input(cpptraj_template, folder, gap)
        outputfile = os.path.join(folder, "complex.pdb")
        print("Fix pdb")
        filename = fix_pdb(folder, filename)
        print(filename + " created")
        assert os.path.exists(filename)

        print("Running tleap:")
        tleap_file = template_tleap(leap_template, filename, outputfile, folder)
        tleap_log = os.path.join(folder, "_tleap.log")
        os.system("tleap -sf {} > {}".format(tleap_file, tleap_log))
        assert os.path.exists(outputfile)
        print("tleap -sf {} > {}".format(tleap_file, tleap_log))

        print("Runningcpptraj:")
        out_file = template_cpptraj(cpptraj_clean, outputfile, folder)
        cpptraj_log = os.path.join(folder, "_cpptraj.log")
        os.system("cpptraj -i {} > {}".format(out_file, cpptraj_log))
        print("cpptraj -i {} > {}".format(out_file, cpptraj_log))

        print("Comparing Fps")
        fp = get_fp(folder)
        comparation, extremes = compare_fp(fp, MD_folder)
        return comparation

def main(pdbs, cpptraj_template="cpptraj_all_fp.txt", leap_template="leap.in",
    MD_folder="fp_BU_MD", gap=1):
    #Define variables
    files = glob.glob(pdbs)

    # Get all input files
    comparations = [extract_info_file(f, cpptraj_template, leap_template, MD_folder, gap) for f in files]

    # Sum number of fps
    score = [sum(comp) for comp in comparations]

    # Sort files by score
    idxs = np.argsort(score)[::-1]
    files_sorted = np.array(files)[idxs]
    scores_sorted = np.array(score)[idxs]

    # Print results
    results = {}
    for file, score in zip(files_sorted, scores_sorted):
        print(file, score)
        results[file]=score

    # Save results
    np.save("score", results) 

    return results
    
        


if __name__ == "__main__":
    #main("*5C*.pdb", gap=+4)
    #main("*6D*.pdb", gap=-1)
    #main("*epoch*.pdb", gap=-1)
    main("*snap1*.pdb", gap=0)
