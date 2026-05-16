# simulate_genotypes.py – Script to simulate genotypes (originally by Iain Mathieson)
#
# Copyright 2018 Iain Mathieson
# Modifications Copyright 2025 Your Name (University of XYZ)
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#
# NOTICE: This file is adapted from the original work by Iain Mathieson (GitHub repository: https://github.com/mathii/gdc).
# Certain functions and logic are reused here with modifications.
# Modified by Bing Li on 2025-04-10 to integrate into this project and to add new features as described:

# The original code is written in python v2 --> changes kept in comments 
# The below code is modified to python3
# In addition, the script extracts the gdc.open2 function from: https://github.com/mathii/gdc/blob/master/gdc.py 

# Below are comments are the original codes: 

# Convert a vcf file to eigenstrat format
# removes multi-alleleic and indel sites. 
# usage: python vcf2eigenstrat.py -v vcf_file.vcf(.gz) -o out_root
# will generate out_root.[snp,ind,geno].
# removed multiallelic sites and indels
# Deals with haploid cases including mixed haploid/diplod like X as well. 
# -i option is a .ind file to get population names and sex. 

from __future__ import division
import sys, getopt, gzip 

################################################################################ 
# The below function is extracted from: https://github.com/mathii/gdc/blob/master/gdc.py 
def open2(file, mode="r"):
	"""
	Open a file, or a gzipped file if it ends in .gz
	"""
	if file[-3:]==".gz":
		return gzip.open(file, mode)
	else:
		return open(file, mode)

################################################################################

def parse_options():
    """
    Options are described by the help() function
    """
    options ={ "vcf":None, "out":"out", "ref":None, "indAsPop":False, "indmap":None  }
	
    try:
        opts, args = getopt.getopt(sys.argv[1:], "v:o:r:i:", ["vcf", "out", "ref", "indmap", "indAsPop"])
        print(opts, args)
        # print opts, args --> original py2
    except Exception as err:
        print (str(err)) # change from print str(err)
        sys.exit()

    for o, a in opts:
        print (o,a) # changed from print o,a
        if o in ["-v","--vcf"]:         options["vcf"] = a
        if o in ["-r","--ref"]:         options["ref"] = a
        if o in ["-i","--ind"]:         options["indmap"] = a
        if o in ["--indAsPop"]:         options["indAsPop"] = True
        elif o in ["-o","--out"]:       options["out"] = a

    # print "found options:"
    # print options
    print("found options:")
    print(options)

    return options

################################################################################

def main(options):
    """
    Convert vcf to eigenstrat format (ind, snp and geno files)
    """
    vcf=open2(options["vcf"])
    snp, ind, geno = [open(options["out"]+x, "w") for x in [".snp", ".ind", ".geno"]]
    removed={"multiallelic":0, "indel":0}
    count=0
    
    if options["indmap"]:
        pop_map={}
        sex_map={}
        ind_map_file=open(options["indmap"], "r")
        for line in ind_map_file:
            bits=line[:-1].split()
            pop_map[bits[0]]=bits[2]
            sex_map[bits[0]]=bits[1]
        ind_map_file.close()
    
    for line in vcf:
        if line[:2]=="##":				  # Comment line
            next
        elif line[:6]=="#CHROM":			  # Header line
            inds=line.split()[9:]
            if options["ref"]:
                ind.write(options["ref"]+"\tU\tREF\n")
            
            if options["indmap"]:
                for indi in inds:
                    ind.write(indi+"\t"+sex_map.get(indi, "U")+"\t"+pop_map.get(indi, "POP")+"\n")
            elif options["indAsPop"]:
                for indi in inds:
                    ind.write(indi+"\tU\t"+indi+"\n")
            else:
                for indi in inds:
                    ind.write(indi+"\tU\tPOP\n")
                   
        else:							  # data
            bits=line.split()
            
            #----------Changed by Bing Li 2025-11-16: Handle phylogenetic VCF with * as first ALT----------#
            # Handle VCF with * as first ALT (phylogenetic VCF format)
            # In VCF: GT indices refer to alleles: 0=REF, 1=first ALT, 2=second ALT, etc.
            # * represents missing/ancestral, real variants start at index 2 in ALT field
            alt_alleles = bits[4].split(",")
            
            # Check if this is a phylogenetic VCF with * as first ALT
            if alt_alleles[0] == "*":
                # Filter to real ALT alleles (skip the *) - only SNPs
                real_alts = [a for a in alt_alleles[1:] if len(a) == 1 and len(bits[3]) == 1]
                
                if len(real_alts) == 0:
                    # No valid biallelic SNP variants, skip this site
                    continue
                
                # For multiallelic sites, we'll record presence/absence of ANY derived allele
                # Use only the first real ALT for the SNP file
                first_real_alt = real_alts[0]
                
                if bits[2]==".":
                    bits[2]=bits[0]+":"+bits[1]
                snp.write("    ".join([bits[2], bits[0], "0.0", bits[1], bits[3], first_real_alt])+"\n")
                
                geno_string=""
                if options["ref"]:
                    geno_string="2"
                
                # Decode genotypes for EIGENSTRAT:
                # EIGENSTRAT: 2=homozygous ancestral/ref, 1=heterozygous, 0=homozygous derived, 9=missing
                # VCF GT: 0=REF (ancestral), 1=* (missing), 2+=variant alleles (derived)
                for gt_field in bits[9:]:
                    gt = gt_field.split(":")[0]
                    if gt == "0":
                        geno_string += "2"  # REF is ancestral (homozygous ancestral)
                    elif gt == "1":
                        geno_string += "9"  # * is missing
                    elif gt == "2":
                        geno_string += "0"  # First real ALT is derived (homozygous derived)
                    elif gt.isdigit() and int(gt) > 2:
                        # Other ALT alleles - treat as derived if only one real ALT, else missing
                        if len(real_alts) == 1:
                            geno_string += "9"  # Wrong allele, treat as missing
                        else:
                            geno_string += "0"  # Another derived allele, lump together
                    else:
                        geno_string += "9"  # Missing or invalid
                
                geno.write(geno_string+"\n")
                count += 1
                #----------End of change by Bing Li 2025-11-16----------#
            else:
                # Standard biallelic VCF
                if "," in bits[4] or len(bits[3]) != 1 or len(bits[4]) != 1:
                    # Skip multiallelic or indels
                    removed["multiallelic"] += 1
                    continue
                
                if bits[2]==".":
                    bits[2]=bits[0]+":"+bits[1]
                snp.write("    ".join([bits[2], bits[0], "0.0", bits[1], bits[3], bits[4]])+"\n")
                
                geno_string=""
                if options["ref"]:
                    geno_string="2"
                
                for gt in bits[9:]:
                    geno_string += decode_gt_string(gt)
                
                geno.write(geno_string+"\n")
                count += 1 

    [f.close for f in [ind, snp, geno]]

    # print "Done. Wrote "+str(count) + " sites"
    # print "Excluded " + str(sum(removed.values())) + " sites"
    # for key in removed:
    #     print "Excluded " + str(removed[key]) + " " + key

    
    print ("Done. Wrote "+str(count) + " sites")
    print ("Excluded " + str(sum(removed.values())) + " sites")
    for key in removed:
        print ("Excluded " + str(removed[key]) + " " + key) 

    return

################################################################################

def decode_gt_string(gt_string):
    """
    Tries to work out the genotype from a vcf genotype entry. 9 for missing [or not in {0,1,2}]
    """
    gt=gt_string.split(":")[0]
    if len(gt)==1:
        if gt=="0":            
            return "2"
            # return "1"
        elif gt=="1":
            return "0"
        else:
            return "9"
    elif len(gt)==3:
        if gt[0]=="0" and gt[2]=="0":
            return "2"
        if gt[0]=="0" and gt[2]=="1":
            return "1"
        if gt[0]=="1" and gt[2]=="0":
            return "1"
        if gt[0]=="1" and gt[2]=="1":
            return "0"
        else:
            return "9"

    raise Exception("Unknown genotype: "+gt)

################################################################################

if __name__=="__main__":
	options=parse_options()
	main(options)
	
