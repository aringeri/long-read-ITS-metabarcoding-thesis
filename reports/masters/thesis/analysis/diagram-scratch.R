DiagrammeR::grViz("digraph {
  graph [layout = dot, rankdir = LR]

  node [shape = rectangle, peripheries=2]
  rec1 [label = 'Adapter trimming\n(dorado)']
  rec2 [label = 'Primer trimming\n(cutadapt)']
  rec3 [label = 'Barcode region extraction\n(ITSxpress)']
  rec4 [label = 'Quality filtering\n(chopper)']
  rec5 [label = 'Chimera filtering\n(vsearch)']

  node [shape = rectangle, peripheries=1]
  rec7 [label = 'UNITE DB']

  node [shape = plaintext, peripheries=0]
  rec0 [label = 'raw demultiplexed reads\nfastq']
  rec6 [label = 'Quality filtered\nfull ITS sequences']


  # edge definitions with the node IDs
  rec0 -> rec1 -> rec2 -> rec3
  rec3 -> rec4 [ label = 'Full ITS']
  rec7 -> rec5
  rec4 -> rec5 -> rec6
  }",
  height = 200)


DiagrammeR::grViz("digraph {
  graph [layout = dot, rankdir = LR]

  node [shape = rectangle, peripheries=1]

  vsearch [ label = '97% identity clustering\n(vsearch)' ]
  kmer [ label = 'k-mer frequencies' ]
  umap [ label = 'Dimension reduction\n(UMAP)' ]
  hdbscan [ label = 'Clustering\n(HDBSCAN)' ]
  mostA [ label = 'Most abundant sequence\n(vsearch)' ]

  dnabarcoder [ label = 'Taxonomic Assignment\n(dnabarcoder)' ]
  unite [ label = 'UNITE DB' ]

  node [shape = rectangle, peripheries=2]
  subs [label = 'Scenario 1\nEven abundance\n(6 libs, 5 reps)']

  node [shape = plaintext, peripheries=0]
  qc [label = 'Quality filtered\nfull ITS sequences']


  # edge definitions with the node IDs
  qc -> subs
  subs -> vsearch
  subs -> kmer -> umap -> hdbscan -> mostA

  vsearch -> dnabarcoder [ label = 'centroids' ]
  mostA -> dnabarcoder [ label = 'most abundant' ]
  unite -> dnabarcoder

  }",
  height = 200)