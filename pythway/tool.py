from . import _CURRENT_WORKING_DIRECTORY, pd, os, datetime
from . import requests as req


def get_kegg_pathway_geneset(org:str=None, force_update:bool=None):
    """
    使用kegg的restAPI读取通路ID及其对应基因ID
    :param org: 物种名称：小鼠mmu，人hsa
    :param force_update: 强制更新缓存文件，默认为真
    :return: pd.Dataframe
    """
    if org is None or force_update is None:
        raise ValueError("Please check setup of org and force_update")
    url = f"https://rest.kegg.jp/link/{org}/pathway"
    # 保存爬取的pathway和gene_ID对应结果
    path = f"{_CURRENT_WORKING_DIRECTORY}/db/{org}_pathway_gene_id.tsv"
    # 如果要求强制更新或者文件不存在
    if not os.path.exists(path) or force_update:
        with open(path, "w", encoding="utf-8") as f:
            # 开始爬取
            pathway2gene_id = req.get(url).content.decode("utf-8").replace("path:","")
            print(f"Trying to get \"{url}\"")
            # 共两列，从左至右依次为pathway_id，gene_id
            f.write("pathway_id\tgene_id\n" + pathway2gene_id)

    return pd.read_csv(path, sep="\t")

def get_kegg_gene_info(org:str=None, force_update:bool=None):
    """
    使用kegg的restAPI读取基因ID及其对应信息
    :param org: 物种名称：小鼠mmu，人hsa
    :param force_update: 强制更新缓存文件，默认为真
    :return: pd.Dataframe
    """
    if org is None or force_update is None:
        raise ValueError("Please check setup of org and force_update")
    # 将gene id转换为gene symbol
    url = f"https://rest.kegg.jp/list/{org}"
    # 保存爬取的gene list
    path = f"{_CURRENT_WORKING_DIRECTORY}/db/{org}_gene_list.tsv"
    # 如果要求强制更新或者文件不存在
    if not os.path.exists(path) or force_update:
        with open(path, "w+", encoding="utf-8") as f:
            # 开始爬取
            gene_list = req.get(url).content.decode("utf-8")
            print(f"Trying to get \"{url}\"")
            # 共四列，从左至右依次为gene_id，gene_type，position，description
            f.write("gene_id\tgene_type\tposition\tdescription\n" + gene_list)
        # 提取基因名
        gene_info = pd.read_csv(path, sep="\t")
        gene_info.loc[:, "gene_symbol"] = gene_info.loc[:, "description"].str.extract(r"^([A-Za-z\.\d-]+)(?=[;,])").iloc[:,0]
        gene_info.to_csv(path, sep="\t", encoding="utf-8", index=False)

    return pd.read_csv(path, sep="\t")

def get_kegg_pathway_description(org:str=None, force_update:bool=None):
    """
    使用kegg的restAPI读取通路ID及其对应信息
    :param org: 物种名称：小鼠mmu，人hsa
    :param force_update: 强制更新缓存文件，默认为真
    :return: pd.Dataframe
    """
    if org is None or force_update is None:
        raise ValueError("Please check setup of org and force_update")
    url = f"https://rest.kegg.jp/list/pathway/{org}"
    # 保存爬取的pathway和对应描述
    path = f"{_CURRENT_WORKING_DIRECTORY}/db/{org}_pathway_description.tsv"
    # 如果要求强制更新或者文件不存在
    if not os.path.exists(path) or force_update:
        with open(path, "w", encoding="utf-8") as f:
            # 开始爬取
            pathway_description = req.get(url).content.decode("utf-8")
            print(f"Trying to get \"{url}\"")
            # 共两列，从左至右依次为pathway_id，gene_id
            f.write("pathway_id\tdescription\n" + pathway_description)

    return pd.read_csv(path, sep="\t")


def get_kegg_pathway_gmt(org:str=None, force_update:bool=None, save_path:str=None):
    """
    生成GSEA使用的gmt文件
    :param org: 物种名称：小鼠mmu，人hsa
    :param force_update: 强制更新缓存文件，默认为真
    :param save_path: gmt文件保存路径
    :return: 无输出
    """
    pathway_gene_id_set = get_kegg_pathway_geneset(org=org, force_update=force_update)
    gene_info = get_kegg_gene_info(org=org, force_update=force_update)
    pathway_description = get_kegg_pathway_description(org=org, force_update=force_update)
    # 匹配基因名
    temp = pd.merge(pathway_gene_id_set, gene_info, on="gene_id", how="left")
    # 如果基因名有NA值，用gene_id
    temp["gene_symbol"] = temp["gene_symbol"].fillna(temp["gene_id"])

    # 连接基因名
    gmt = temp.loc[:, ["pathway_id", "gene_symbol"]].groupby('pathway_id')["gene_symbol"].apply(lambda x: "\t".join(x))
    # 将基因名与通路描述合并
    gmt = pd.merge(pathway_description, gmt, on="pathway_id", how="left")
    # 检查保存路径
    if save_path is None:
        save_path = f"{_CURRENT_WORKING_DIRECTORY}/db/"
    # 保存gmt
    with open(os.path.join(save_path, f"{org}_kegg_v{datetime.datetime.now().strftime("%Y.%m")}_pathway.gmt"),
              "w", encoding="utf-8") as f:
        for _,row in gmt.iterrows():
            f.write(f"{row.iloc[0]}\t{row.iloc[1]}\t{row.iloc[2]}\n")

    return 0
