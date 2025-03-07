import pythway as pw

if __name__ == "__main__":
    pathway_gene_id_set = pw.get_kegg_pathway_geneset(org="mmu", force_update=False)
    gene_info = pw.get_kegg_gene_info(org="mmu", force_update=False)
    pathway_description = pw.get_kegg_pathway_description(org="mmu", force_update=False)
    pw.get_kegg_pathway_gmt(org="mmu", force_update=False, save_path=r"D:\Python\bio_informatics\SAH20240912\data")
    print("end")