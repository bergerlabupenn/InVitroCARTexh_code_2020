#!/usr/bin/python
# Compares KLRB1 expression of TCR-fingerprinted CD8+ cells at day 0 and day 28 CS
# Uses cell barcodes to count cells
# This is necessary because two cells can have different TCRs
# Expects filtered coverage annotation (FCA) files in a folder called "FCAs"
# Expects the CSV files with gene expression data in a folder called
# "Outs/PATIENT.Day0" or "Outs/PATIENT.Day28CS"
# Greg Donahue, 01-30-2020

import sys, string

klrb1 = "ENSG00000111796"
cd8a = "ENSG00000153563"

def main(args):
    
    # For each sample
    for patient in [ "ND150", "ND538"]:
        
        # Load common TCRs and dictionaries mapping TCR->[ cell barcodes ]
        day0_tcrs = loadFCA("FCAs/"+patient+".Day0.FCA.csv")
        day28cs_tcrs = loadFCA("FCAs/"+patient+".Day28cs.FCA.csv")
        common_tcrs = list(set(list(day0_tcrs.keys())).intersection(set(list(day28cs_tcrs.keys()))))
        
        # Load dictionaries mapping cell barcode->[ genes expressed ]
        day0_rna = loadRNA("Outs/"+patient+".Day0/filtered_feature_bc_matrix/"+patient+".Day0.csv")
        day28cs_rna = loadRNA("Outs/"+patient+".Day28CS/filtered_feature_bc_matrix/"+patient+".Day28CS.csv")

        # Load all cell barcodes for cells expressing crucical KLRB1 and CD8A markers
        day0_expressing_klrb1 = collectCellsExpressing(day0_rna, klrb1)
        day28cs_expressing_klrb1 = collectCellsExpressing(day28cs_rna, klrb1)
        day0_expressing_cd8a = collectCellsExpressing(day0_rna, cd8a)
        day28cs_expressing_cd8a = collectCellsExpressing(day28cs_rna, cd8a)

        # Make a dictionary going from day 0 to TCR (just for common TCRs)
        day0_to_tcr = reverseDictionary(day0_tcrs, common_tcrs)

        # Identify persistent, CD8+ cells which gain KLRB+ in day 28 CS
        # Strategy: start with day 0 cell barcodes for cells with common TCRs and CD8A expression, then
        #      find the TCRs, and ask if any cell at day 28 CS attached to those TCRs has CD8A expression and/or
        #      KLRB1 expression
        #      We are negotiating the issue of doublets by guessing that when one cell exhibits two TCRs, it's really
        #      two cells
        constitutive_on, constitutive_off, turned_on, turned_off = list(), list(), list(), list()
        day0_cd8_persistent_cells = list(set(list(day0_to_tcr.keys())).intersection(set(day0_expressing_cd8a)))
        partners = dict()
        for C in day0_cd8_persistent_cells:
            day28_cs_partners = list()
            for TCR in day0_to_tcr[C]: day28_cs_partners.extend(day28cs_tcrs[TCR])
            day28_cs_partners = list(set(day28_cs_partners).intersection(set(day28cs_expressing_cd8a)))
            if len(day28_cs_partners) == 0: continue
            partners[C] = day28_cs_partners
            if C in day0_expressing_klrb1:
                if len(set(day28_cs_partners).intersection(set(day28cs_expressing_klrb1))) > 0:
                    constitutive_on.append(C)
                else: turned_off.append(C)
            else:
                if len(set(day28_cs_partners).intersection(set(day28cs_expressing_klrb1))) > 0:
                    turned_on.append(C)
                else: constitutive_off.append(C)
        print("Report for", patient+":")
        print("Constitutive On:", str(len(constitutive_on)))
        for C in constitutive_on:
            print("\t"+C+"\t"+",".join(day0_to_tcr[C])+"\t"+",".join(partners[C]))
        print("Turned On:", str(len(turned_on)))
        for C in turned_on:
            print("\t"+C+"\t"+",".join(day0_to_tcr[C])+"\t"+",".join(partners[C]))
        print("Turned Off:", str(len(turned_off)))
        for C in turned_off:
            print("\t"+C+"\t"+",".join(day0_to_tcr[C])+"\t"+",".join(partners[C]))
        print("Constitutive Off:", str(len(constitutive_off)))
        for C in constitutive_off:
            print("\t"+C+"\t"+",".join(day0_to_tcr[C])+"\t"+",".join(partners[C]))

def reverseDictionary(data, filter_list):
    r_dict = dict()
    for tcr in filter_list:
        for cell in data[tcr]:
            try: r_dict[cell].append(tcr)
            except Exception as e: r_dict[cell] = [ tcr ]
    return r_dict

def collectCellsExpressing(data, gene):
    r_cells = list()
    for cell in list(data.keys()):
        if len(set([ gene ]).intersection(set(data[cell]))) > 0: r_cells.append(cell)
    return r_cells

def loadRNA(filename):
    r_cell_data = dict()
    with open(filename) as f:
        cells = f.readline().rstrip().split(",")[1:]
        for C in cells: r_cell_data[C] = list()
        for line in f:
            t = line.rstrip().split(",")
            for i in list(range(len(t)-1)):
                if int(t[i+1]) > 0: r_cell_data[cells[i]].append(t[0])
    return r_cell_data

def loadFCA(filename):
    r_tcrs = dict()
    with open(filename) as f:
        next(f)
        for line in f:
            t = line.rstrip().split(",")
            if t[10] != "True": continue
            if t[12] == "None": continue
            try: r_tcrs[t[12]].append(t[0])
            except Exception as e: r_tcrs[t[12]] = [ t[0] ]
    return r_tcrs
        

if __name__ == "__main__": main(sys.argv)
