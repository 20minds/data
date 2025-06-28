if not os.path.exists('pqtl_ukbb.parquet'):
    urllib.request.urlretrieve('https://github.com/20minds/data/raw/refs/heads/main/bio/pqtl_ukbb.parquet', './pqtl_ukbb.parquet')
df = pd.read_parquet('pqtl_ukbb.parquet')
df.loc[0, :].to_dict()
{'rs_id': 'rs41265213', 'gene_id': 'ENSG00000160712', 'cell_type_name': 'SCALLOP_2020-UBERON_0001969', 'qtl_score': 5.13715334001706, 'gene_name': 'IL6R'}

if not os.path.exists('eqtl_ukbb.parquet'):
    urllib.request.urlretrieve('https://github.com/20minds/data/raw/refs/heads/main/bio/eqtl_ukbb.parquet', './eqtl_ukbb.parquet')
df = pd.read_parquet('eqtl_ukbb.parquet')
df.loc[0, :].to_dict()
{'rs_id': 'rs111972115', 'gene_id': 'ENSG00000134250', 'cell_type_name': 'CEDAR-MONOCYTE_CD14', 'qtl_score': 4.851273983601853, 'gene_name': 'NOTCH2'}

if not os.path.exists('gwas_catalog.parquet'):
    urllib.request.urlretrieve('https://github.com/20minds/data/raw/refs/heads/main/bio/gwas_catalog.parquet', './gwas_catalog.parquet')
df = pd.read_parquet('gwas_catalog.parquet')
df.loc[0, :].to_dict()
{'DATE ADDED TO CATALOG': '2018-07-30', 'PUBMEDID': 27618447, 'FIRST AUTHOR': 'Surendran P', 'DATE': '2016-10-01', 'JOURNAL': 'Nat Genet', 'LINK': 'www.ncbi.nlm.nih.gov/pubmed/27618447', 'STUDY': 'Trans-ancestry meta-analyses identify rare and common variants associated with blood pressure and hypertension.', 'DISEASE/TRAIT': 'systolic blood pressure', 'INITIAL SAMPLE SIZE': 'up to 165,276 European ancestry individuals, up to 27,487 South Asian ancestry individuals', 'REPLICATION SAMPLE SIZE': 'up to 125,713 European ancestry individuals, up to 2,641 South Asian ancestry individuals, 4,632 Hispanic individuals, 22,077 African American individuals', 'REGION': '3q26.2', 'CHR_ID': '3', 'CHR_POS': '169383111', 'REPORTED GENE(S)': 'MECOM', 'MAPPED_GENE': 'MECOM', 'UPSTREAM_GENE_ID': None, 'DOWNSTREAM_GENE_ID': None, 'SNP_GENE_IDS': 'ENSG00000085276', 'UPSTREAM_GENE_DISTANCE': nan, 'DOWNSTREAM_GENE_DISTANCE': nan, 'STRONGEST SNP-RISK ALLELE': 'rs448378-A', 'SNPS': 'rs448378', 'MERGED': 0, 'SNP_ID_CURRENT': '448378.0', 'CONTEXT': 'intron_variant', 'INTERGENIC': 0.0, 'RISK ALLELE FREQUENCY': '0.5261', 'P-VALUE': 8e-09, 'PVALUE_MLOG': 8.096910013008056, 'P-VALUE (TEXT)': None, 'OR or BETA': 0.0196, '95% CI (TEXT)': 'mmHg decrease', 'PLATFORM [SNPS PASSING QC]': 'Illumina [242296]', 'CNV': 'N'}

if not os.path.exists('proximity_label_ms.parquet'):
    urllib.request.urlretrieve('https://github.com/20minds/data/raw/refs/heads/main/bio/proximity_label_ms.parquet', './proximity_label_ms.parquet')
df = pd.read_parquet('proximity_label_ms.parquet')
df.loc[0, :].to_dict()
{'interaction_id': 913336, 'gene_a_id': 'ENSG00000206560', 'gene_b_id': 'ENSG00000103194', 'experimental_system_type': 'physical', 'pubmed_id': 'PUBMED:24255178', 'organism_id_a': 9606, 'organism_id_b': 9606, 'throughput_type': 'High Throughput', 'experimental_score': 1.0}

if not os.path.exists('co_fractionation.parquet'):
    urllib.request.urlretrieve('https://github.com/20minds/data/raw/refs/heads/main/bio/co_fractionation.parquet', './co_fractionation.parquet')
df = pd.read_parquet('co_fractionation.parquet')
df.loc[0, :].to_dict()
{'interaction_id': 740709, 'gene_a_id': 'ENSG00000160201', 'gene_b_id': 'ENSG00000063244', 'experimental_system_type': 'physical', 'pubmed_id': 'PUBMED:22939629', 'organism_id_a': 9606, 'organism_id_b': 9606, 'throughput_type': 'High Throughput', 'experimental_score': 0.874}

if not os.path.exists('reconstituted_complex.parquet'):
    urllib.request.urlretrieve('https://github.com/20minds/data/raw/refs/heads/main/bio/reconstituted_complex.parquet', './reconstituted_complex.parquet')
df = pd.read_parquet('reconstituted_complex.parquet')
df.loc[0, :].to_dict()
{'interaction_id': 2841862, 'gene_a_id': 'ENSG00000169083', 'gene_b_id': 'ENSG00000097007', 'experimental_system_type': 'physical', 'pubmed_id': 'PUBMED:24728074', 'organism_id_a': 9606, 'organism_id_b': 9606, 'throughput_type': 'High Throughput', 'experimental_score': 2.13}

if not os.path.exists('affinity_capture_rna.parquet'):
    urllib.request.urlretrieve('https://github.com/20minds/data/raw/refs/heads/main/bio/affinity_capture_rna.parquet', './affinity_capture_rna.parquet')
df = pd.read_parquet('affinity_capture_rna.parquet')
df.loc[0, :].to_dict()
{'interaction_id': 637589, 'gene_a_id': 'YDR515W', 'gene_b_id': 'YAL030W', 'experimental_system_type': 'physical', 'pubmed_id': 'PUBMED:22271760', 'organism_id_a': 559292, 'organism_id_b': 559292, 'throughput_type': 'High Throughput', 'experimental_score': 0.736666667}

if not os.path.exists('genebass_missense_LC_filtered.parquet'):
    urllib.request.urlretrieve('https://github.com/20minds/data/raw/refs/heads/main/bio/genebass_missense_LC_filtered.parquet', './genebass_missense_LC_filtered.parquet')
df = pd.read_parquet('genebass_missense_LC_filtered.parquet')
df.loc[0, :].to_dict()
{'annotation': 'missense_LC', 'Pvalue': 0.85434, 'Pvalue_Burden': 0.69637, 'Pvalue_SKAT': 0.87135, 'BETA_Burden': -0.0046974, 'SE_Burden': 0.012037, 'gene': 'TSPAN6', 'pheno_description': 'Treatment/medication code; co-dydramol'}

if not os.path.exists('genebass_synonymous_filtered.parquet'):
    urllib.request.urlretrieve('https://github.com/20minds/data/raw/refs/heads/main/bio/genebass_synonymous_filtered.parquet', './genebass_synonymous_filtered.parquet')
df = pd.read_parquet('genebass_synonymous_filtered.parquet')
df.loc[0, :].to_dict()
{'annotation': 'synonymous', 'Pvalue': 0.72199, 'Pvalue_Burden': 0.85198, 'Pvalue_SKAT': 0.52215, 'BETA_Burden': 0.0034703, 'SE_Burden': 0.018598, 'gene': 'TSPAN6', 'pheno_description': 'Eye problems/disorders; Macular degeneration'}

if not os.path.exists('genebass_pLoF_filtered.parquet'):
    urllib.request.urlretrieve('https://github.com/20minds/data/raw/refs/heads/main/bio/genebass_pLoF_filtered.parquet', './genebass_pLoF_filtered.parquet')
df = pd.read_parquet('genebass_pLoF_filtered.parquet')
df.loc[0, :].to_dict()
{'annotation': 'pLoF', 'Pvalue': 0.59145, 'Pvalue_Burden': 0.96398, 'Pvalue_SKAT': 0.42053, 'BETA_Burden': -0.00053797, 'SE_Burden': 0.011913, 'gene': 'DPM1', 'pheno_description': 'Mean L1 in fornix on FA skeleton'}

if not os.path.exists('synthetic_growth_defect.parquet'):
    urllib.request.urlretrieve('https://github.com/20minds/data/raw/refs/heads/main/bio/synthetic_growth_defect.parquet', './synthetic_growth_defect.parquet')
df = pd.read_parquet('synthetic_growth_defect.parquet')
df.loc[0, :].to_dict()
{'interaction_id': 299191, 'gene_a_id': 'SPBC30D10.04', 'gene_b_id': 'SPAC14C4.13', 'experimental_system_type': 'genetic', 'pubmed_id': 'PUBMED:18931302', 'organism_id_a': 284812, 'organism_id_b': 284812, 'throughput_type': 'High Throughput', 'experimental_score': -366.0}

if not os.path.exists('genetic_interaction.parquet'):
    urllib.request.urlretrieve('https://github.com/20minds/data/raw/refs/heads/main/bio/genetic_interaction.parquet', './genetic_interaction.parquet')
df = pd.read_parquet('genetic_interaction.parquet')
df.loc[0, :].to_dict()
{'interaction_id': 206363, 'gene_a_id': 'YCR011C', 'gene_b_id': 'YCL025C', 'experimental_system_type': 'genetic', 'pubmed_id': 'PUBMED:16269340', 'organism_id_a': 559292, 'organism_id_b': 559292, 'throughput_type': 'High Throughput', 'experimental_score': -5.6431}

if not os.path.exists('two_hybrid.parquet'):
    urllib.request.urlretrieve('https://github.com/20minds/data/raw/refs/heads/main/bio/two_hybrid.parquet', './two_hybrid.parquet')
df = pd.read_parquet('two_hybrid.parquet')
df.loc[0, :].to_dict()
{'interaction_id': 269293, 'gene_a_id': 'ENSG00000114395', 'gene_b_id': 'ENSG00000109103', 'experimental_system_type': 'physical', 'pubmed_id': 'PUBMED:16169070', 'organism_id_a': 9606, 'organism_id_b': 9606, 'throughput_type': 'High Throughput', 'experimental_score': 1.0}

if not os.path.exists('trait.parquet'):
    urllib.request.urlretrieve('https://github.com/20minds/data/raw/refs/heads/main/bio/trait.parquet', './trait.parquet')
df = pd.read_parquet('trait.parquet')
df.loc[0, :].to_dict()
{'pheno_name': '120_1', 'description': 'Birth weight known; Yes - pounds and ounces'}

if not os.path.exists('affinity_capture_ms.parquet'):
    urllib.request.urlretrieve('https://github.com/20minds/data/raw/refs/heads/main/bio/affinity_capture_ms.parquet', './affinity_capture_ms.parquet')
df = pd.read_parquet('affinity_capture_ms.parquet')
df.loc[0, :].to_dict()
{'interaction_id': 241421, 'gene_a_id': 'ENSG00000169217', 'gene_b_id': 'ENSG00000144028', 'experimental_system_type': 'physical', 'pubmed_id': 'PUBMED:17353931', 'organism_id_a': 9606, 'organism_id_b': 9606, 'throughput_type': 'High Throughput', 'experimental_score': 1.0}

if not os.path.exists('sqtl_ukbb.parquet'):
    urllib.request.urlretrieve('https://github.com/20minds/data/raw/refs/heads/main/bio/sqtl_ukbb.parquet', './sqtl_ukbb.parquet')
df = pd.read_parquet('sqtl_ukbb.parquet')
df.loc[0, :].to_dict()
{'rs_id': 'rs12756687', 'gene_id': 'ENSG00000265491', 'cell_type_name': 'GTEx-sQTL-Brain_Caudate_basal_ganglia', 'qtl_score': 3.19846622732464, 'gene_name': 'RNF115'}

if not os.path.exists('gene_info.parquet'):
    urllib.request.urlretrieve('https://github.com/20minds/data/raw/refs/heads/main/bio/gene_info.parquet', './gene_info.parquet')
df = pd.read_parquet('gene_info.parquet')
df.loc[0, :].to_dict()
{'gene_id': 'ENSG00000228037', 'transcript_id': 'ENST00000424215', 'chr': '1', 'gene_start': 2581560, 'gene_end': 2584533, 'strand': 1, 'transcript_start': 2581560, 'transcript_end': 2584533, 'tss': 2581560, 'transcript_is_canonical': 1.0, 'gene_name': None, 'percentage_gene_gc_content': 51.11, 'gene_type': 'lncRNA'}

if not os.path.exists('gtex_tissue_gene_tpm.parquet'):
    urllib.request.urlretrieve('https://github.com/20minds/data/raw/refs/heads/main/bio/gtex_tissue_gene_tpm.parquet', './gtex_tissue_gene_tpm.parquet')
df = pd.read_parquet('gtex_tissue_gene_tpm.parquet')
df.loc[0, :].to_dict()
{'Description': 'ENSG00000186092', 'Tissue': 'Adipose - Subcutaneous', 'Expression': 0.0453961, 'Gene': 'OR4F5'}

if not os.path.exists('variant_table.parquet'):
    urllib.request.urlretrieve('https://github.com/20minds/data/raw/refs/heads/main/bio/variant_table.parquet', './variant_table.parquet')
df = pd.read_parquet('variant_table.parquet')
df.loc[0, :].to_dict()
{'RS': 'rs116587930', 'ID': '1:727841_G_A', 'CHR': 1, 'POS': 727841, 'A1': 'G', 'A2': 'A', 'MAF': 0.0507035}

if not os.path.exists('dosage_growth_defect.parquet'):
    urllib.request.urlretrieve('https://github.com/20minds/data/raw/refs/heads/main/bio/dosage_growth_defect.parquet', './dosage_growth_defect.parquet')
df = pd.read_parquet('dosage_growth_defect.parquet')
df.loc[0, :].to_dict()
{'interaction_id': 908565, 'gene_a_id': 'YBL009W', 'gene_b_id': 'YPR058W', 'experimental_system_type': 'genetic', 'pubmed_id': 'PUBMED:22282571', 'organism_id_a': 559292, 'organism_id_b': 559292, 'throughput_type': 'High Throughput', 'experimental_score': -0.504}

if not os.path.exists('synthetic_rescue.parquet'):
    urllib.request.urlretrieve('https://github.com/20minds/data/raw/refs/heads/main/bio/synthetic_rescue.parquet', './synthetic_rescue.parquet')
df = pd.read_parquet('synthetic_rescue.parquet')
df.loc[0, :].to_dict()
{'interaction_id': 986367, 'gene_a_id': 'YJR017C', 'gene_b_id': 'YHR114W', 'experimental_system_type': 'genetic', 'pubmed_id': 'PUBMED:24470217', 'organism_id_a': 559292, 'organism_id_b': 559292, 'throughput_type': 'High Throughput', 'experimental_score': -0.2}

if not os.path.exists('synthetic_lethality.parquet'):
    urllib.request.urlretrieve('https://github.com/20minds/data/raw/refs/heads/main/bio/synthetic_lethality.parquet', './synthetic_lethality.parquet')
df = pd.read_parquet('synthetic_lethality.parquet')
df.loc[0, :].to_dict()
{'interaction_id': 818564, 'gene_a_id': 'YLR418C', 'gene_b_id': 'YLR103C', 'experimental_system_type': 'genetic', 'pubmed_id': 'PUBMED:23390603', 'organism_id_a': 559292, 'organism_id_b': 559292, 'throughput_type': 'High Throughput|Low Throughput', 'experimental_score': 4.75e-06}

