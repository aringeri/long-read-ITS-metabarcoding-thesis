DiagrammeR::grViz("digraph {
  graph [layout = dot, label='Quality Filtering of Raw Reads', labelloc='tl', labeljust='l', rankdir = LR, fontsize='30']

  node [shape = rectangle, peripheries=2, fontsize='20']
  rec1 [label = 'Adapter trimming\n(dorado)']
  rec2 [label = 'Primer trimming\n(cutadapt)']
  rec3 [label = 'Barcode region extraction\n(ITSxpress)']
  rec4 [label = 'Quality filtering\n(chopper)']
    subgraph chimera {
    rank = same
    rec5 [label = 'Chimera filtering\n(vsearch)']
    node [shape = rectangle, peripheries=1]
    rec7 [label = 'UNITE DB']
    rec7 -> rec5;
  }

  node [shape = plaintext, peripheries=0]
  rec0 [label = 'Raw demultiplexed reads\nfastq']
  rec6 [label = 'Quality filtered\nfull ITS sequences']


  # edge definitions with the node IDs
  edge [ fontsize='20']
  rec0 -> rec1 -> rec2 -> rec3
  rec3 -> rec4 [ label = 'Full ITS']
  rec4 -> rec5 -> rec6
  }",
  height = 200)


DiagrammeR::grViz("digraph {
  graph [layout = dot, label='Read Clustering', labelloc='tl', labeljust='l', rankdir = LR, fontsize='50']

  node [shape = rectangle, peripheries=1, fontsize='30', width=3.5]

  subgraph {
    subgraph A {
      derep [ label = 'Read dereplication\n(vsearch)' ]
      vsearch [ label = '97% identity clustering\n(vsearch)' ]
    }
    subgraph B {
      kmer [ label = 'k-mer frequencies\n(python)' ]
      umap [ label = 'Dimension reduction\n(UMAP)' ]
      hdbscan [ label = 'Clustering\n(HDBSCAN)' ]
    }
  }
  subgraph {
    rank = same
    mostA [ label = 'Most abundant sequence\n(vsearch)' ]
    consA [ label = 'Consensus sequence\n(polishing)' ]
    centr [ label = 'centroids\n(vsearch)' ]

  }

  node [shape = rectangle, peripheries=2]
  subs [label = 'Scenario 1\nEven abundance\n(6 libs, 5 reps)']

  node [shape = plaintext, peripheries=0]
  qc [label = 'Quality filtered\nfull ITS sequences']

  reps [label = 'Representative sequences\nfor classification' ]

  # edge definitions with the node IDs
  edge [ fontsize = '20']
  qc -> subs
  subs -> derep -> vsearch
  subs -> kmer -> umap -> hdbscan -> mostA
  hdbscan -> consA
  vsearch -> consA

  vsearch -> centr -> reps
  consA -> reps
  mostA -> reps
  #vsearch -> dnabarcoder [ label = 'centroids' ]
  #mostA -> dnabarcoder [ label = 'most abundant' ]
  #unite -> dnabarcoder
  }",
  height = 200)

DiagrammeR::grViz("digraph {
  graph [layout = dot, label='Consensus Sequence Polishing', labelloc='tl', labeljust='l', rankdir = LR, fontsize='50']

  node [shape = rectangle, peripheries=1, fontsize='30', width=3.5]

  subset [ label = 'Subset reads for polishing\n(n=200)' ]
  subgraph {
  fastANI [ label = 'Select draft by max ANI\n(fastANI)' ]
  racon [  label = 'Polish draft\n(racon)' ]
  medaka [ label = 'Polish draft\n(medaka)' ]
  }

  node [shape = plaintext, peripheries=0]
  reads [label = 'Set of reads in same cluster']
  reps [label = 'Representative sequences\nfor classification' ]

  # edge definitions with the node IDs
  edge [ fontsize = '20']
  reads -> subset -> fastANI -> racon -> medaka -> reps [ weight = 10 ]
  subset -> racon:nw
  subset -> medaka:nw
  }",
  height = 200)