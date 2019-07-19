# protein-group-function-prediction
Prediction of protein group function by iterative classification on functional relevance graph

I. Data

10 protein groups involved in Rheumatoid Arthritis (RA).

Data/Nets/RA/
1. Allograft rejection
2. Apoptosis
3. Pathways in cancer
4. Chemokine signaling
5. Jak-STAT
6. Leukocyte migration
7. MAPK signaling
8. Neurotrophin signaling
9. T cell receptor signaling
10. Toll-like rec. signaling

<geneGroupFileName> the uniprot ID of the disease genes/any group of genes
format:: <genename>	<UAC>

<networkName> graph edge weights for a protein * protein graph

II. The GFP directory contains R and perl codes

perl run_full_pipeleine.pl <geneGroupFileName> <networkName>



