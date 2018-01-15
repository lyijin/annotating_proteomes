====================
Annotating proteomes
====================

Unlike most tutorials I write, this one assumes that the reader is familiar with the UNIX shell.

As a sidenote, I've used this pipeline to annotate the gene models of:

1. Aiptasia pallida (Baumgarten et al., 2015)
2. Symbiodinium microadriaticum (Aranda et al., 2016)
3. Amplexidiscus fenestrafer and Discosoma sp. (Wang et al., 2017)
4. Stylophora pistillata (Voolstra et al., 2017)
5. Platygyra daedalea (unpublished)

... which means that this README would hopefully be more informative than the terse sentences one usually writes in the Materials & Methods of scientific papers.

Pipeline summary
----------------
1. Install ``ncbi-blast+``.
2. Download and prepare the holy trinity of protein database sequences: Swiss-Prot, TrEMBL, and nr.
3. Perform three BLASTP searches, then wait for a few weeks.
4. While BLASTP is running, prepare GO term-related files.
5. Once BLASTP is done, parse the BLASTP XMLs.
6. Compile an overall table.

Have you heard of our lord and saviour, UniProt
-----------------------------------------------
Why Swiss-Prot and TrEMBL, you ask? Proteins in these databases are assigned very descriptive annotations--most importantly of which, GO terms. I wanted GO terms assigned to the gene model of the newly sequenced organisms because I could subsequently carry out functional enrichment analyses on lists of genes of interest.

While the major pro of using Swiss-Prot and TrEMBL is the annotations; the major con is that these databases are nowhere near as exhaustive as nr. This is where nr comes in: if a particular gene model has no hits against Swiss-Prot or TrEMBL, it could still have a hit against nr (albeit without GO term annotations, sadly).

The annotation priority for a new gene model would be:

1. If the gene model has a hit against Swiss-Prot, take it. Otherwise,
2. If the gene model has a hit against TrEMBL, take it. Otherwise,
3. If the gene model has a hit against nr, take it. Otherwise, gene model is unannotated.

OK, enough proselytisation: let's get to work.

Installing ncbi-blast+
----------------------
``apt install ncbi-blast+`` for Debian-likes. I don't use ``yum``, but I suspect a similarly-named repo is available.

Test that it works by doing ``blastp -version``. If it's really old, it might be that you're using a really old distro. Get the newest blast+ executables here: ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/

If you did manually download ncbi-blast+, tweak your ``$PATH`` to include the ``bin`` of blast+.

Preparing the holy trinity
--------------------------
Note that these files do get updated over time. In my projects, I have made it a habit to note down when I downloaded these files.

- SwissProt: ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz
- TrEMBL: ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_trembl.fasta.gz

``gunzip`` these files. I like to rename the uncompressed fasta files to ``sprot`` and ``trembl`` respectively (yes, without a ``.fa`` at the end. You'll see why later).

Convert these fasta files into BLAST databases.

``makeblastdb -in sprot -dbtype prot -title "UniProt/Swiss-Prot (Jan 2018)" -parse_seqids``

``makeblastdb -in trembl -dbtype prot -title "UniProt/TrEMBL (Jan 2018)" -parse_seqids``

For nr, this is a bit hacky.

1. Go to ftp://ftp.ncbi.nlm.nih.gov/blast/db/ and look for files in the pattern of nr.*.tar.gz.
2. Make a note of what the largest number is (79, as of Jan 2018. It's insane how fast these numbers have incremented!).
3. Edit ``ftps_get_nr.sh`` line 3, and put in that largest number (``for a in {00..79}``).
4. Run it, you might need to wait for a day or two. Each file is ~500 MB, so all files would be ~40 GB. I hope you've no data caps.
5. Decompress everything (``for a in *.tar.gz; do tar zxvf ${a}; done``). The compressed 40 GB will become 150 GB or so. Make sure you've got loads of free HDD space.
6. Optional cleanup: to free space, feel free to delete the downloaded ``*.tar.gz`` files, ``sprot``, and ``trembl``.

The good thing about this part is that you can do it once every year-ish, and annotate multiple genomes with these files. If you work on a new genome every two years, then... please re-download these files for the newer project.

Running BLASTP searches
-----------------------
Some optional parameters you'll encounter later are not optional, because my downstream scripts require the output files to be XML. It's a choice I made years ago, as ``blastp`` at that time produced a tabular output that truncated annotations, which irritated me to no end. The tabular format is fixed now, but I'm lazy to adapt my scripts. XML it is!

Let's use *Stylophora pistillata* as an example. Its 4-letter code is "spis".

Run these searches in parallel.

``blastp -db sprot -query Spis.genome.annotation.pep.longest.fa -outfmt 5 -max_target_seqs 20 -num_threads 10 -out spis_vs_sprot.blastp.xml``

``blastp -db trembl -query Spis.genome.annotation.pep.longest.fa -outfmt 5 -max_target_seqs 20 -num_threads 10 -out spis_vs_trembl.blastp.xml``

``blastp -db nr -query Spis.genome.annotation.pep.longest.fa -outfmt 5 -max_target_seqs 20 -num_threads 30 -out spis_vs_nr.blastp.xml``

Feel free to tweak ``-query``, ``-num_threads`` and ``-out``, don't touch the other two parameters.

``-outfmt 5`` produces XML files.

``-max_target_seqs 20`` produces at most 20 hits. Why 20, instead of 1 (i.e. best hit only)? It's because in the process of annotating genomes, I realised that some of the best hits... do not have GO terms assigned to them! Rare, but it does happen! Imagine my frustration when I discovered this fact, which led to the writing of ``get_top_hit_with_amigo_annot.py`` (this'll be explained later, don't worry).

So yes, I presented a half-truth earlier for ease of understanding. The priority is *actually*:

1. If the gene model has hits against Swiss-Prot, grab the top 20 hits. Check whether the best hit has GO terms. If not, check second best. If not, check third best... Otherwise, if none of the hits have GO terms assigned to them,
2. If the gene model has hits against TrEMBL, grab the top 20 hits. Check whether the best hit has GO terms. If not, check second best hit. If not, check third best... Otherwise, if none of the hits have GO terms assigned to them,
3. If the gene model has a hit against nr, take it. Otherwise, gene model is unannotated.

From experience, ``blastp`` against nr is much much slower than the other two (hence why I run it with more threads). Do the next section while the ``blastp`` searches are running.

Preparing GO files
------------------
These files contain the mapping of UniProt ID --> GO term. There's a lot of them... which is why you need to download a few more multi-GB files. Use ``wget``.

GO annotation file: http://www.geneontology.org/gene-associations/goa_uniprot_all.gaf.gz

GO term hierarchy: http://purl.obolibrary.org/obo/go/go-basic.obo

Run the shell script ``parse_gp_assoc.sh`` to produce ``goa_uniprot_all.parsed.gaf`` and ``goa_uniprot_all.unique_ids.txt``.

While we're at it, modify line 35 of ``parse_go_obo.py`` (``go_term_file = open('...')``) to point to where you kept ``go-basic.obo``. This script is required later.

Wait for the ``blastp`` searches to finish.

Parsing the XML outputs
-----------------------
When the ``blastp`` searches finish, you should have three files.

- ``spis_vs_sprot.blastp.xml``
- ``spis_vs_trembl.blastp.xml``
- ``spis_vs_nr.blastp.xml``

In the same folder, run

``parse_blast_xml.py --table -e 1e-5 -t 20 spis_vs_sprot.blastp.xml > spis_vs_sprot.t20.tsv``

``parse_blast_xml.py --table -e 1e-5 -t 20 spis_vs_trembl.blastp.xml > spis_vs_trembl.t20.tsv``

``parse_blast_xml.py --table -e 1e-5 -t 1 spis_vs_nr.blastp.xml > spis_vs_nr.t1.tsv``

The ``-e`` flag controls the e-value cutoff. I tend to use 10\ :sup:`-5` (hence ``-e 1e-5``).

These commands parses the XML files into tabulated BLAST results. At this point, you can save space by compressing the XML files (try not to delete them, these files took weeks to produce--they're actually stuff I archive in case things go wrong downstream).

To get the best hit with GO terms, run

``get_top_hit_with_amigo_annot.py spis_vs_sprot.t20.tsv > spis_vs_sprot.tGO.tsv``

``get_top_hit_with_amigo_annot.py spis_vs_trembl.t20.tsv > spis_vs_trembl.tGO.tsv``

The genes in ``*.tGO.tsv`` are guaranteed to have GO terms annotated to it, which circumvents the issue I noticed years ago (of best hits occasionally not having annotated GO terms).

Assigning GO terms
------------------
The two files you need for this step are:

- ``spis_vs_sprot.tGO.tsv``
- ``spis_vs_trembl.tGO.tsv``

Check that you've amended line 35 of ``parse_go_obo.py`` to point to the location of your ``go_basic.obo`` file.

Remember that my example uses "spis" as the 4-letter species code, and it's in all of my files (``spis_*``).

``create_go_annots_sprot_trembl.py spis``

It should produce four files, the most important being ``spis_go_annots.all.tsv``. Files with ``bp``, ``cc`` and ``mf`` nestled within the filenames correspond to genes annotated with "biological process"-, "cellular component"- and "molecular function"-related terms respectively. I personally don't really find them that useful, but if you're interested in only analysing a particular class of GO terms, you could use these files.

Sidenote: the file in the form of ``*_go_annots.all.tsv`` is what I use in my functional enrichment pipeline.

Compiling the overall annotation table
--------------------------------------
Five files are needed here:

- ``spis_vs_sprot.tGO.tsv``
- ``spis_vs_trembl.tGO.tsv``
- ``spis_vs_nr.t1.tsv``
- ``spis_go_annots.all.tsv``
- A FASTA file of your protein gene models. In my project, it was called ``Spis.genome.annotation.pep.longest.fa``

To compile the table, run

``create_top_hit_sprot_trembl_nr.py spis Spis.genome.annotation.pep.longest.fa > spis_tabulated_annots.tsv``

That's it! I usually format the ``.tsv`` file (tab-separated text file) into an Excel sheet, and delete a few unnecessary columns for publication purposes; for normal bioinformatics-y work, I use the ``.tsv`` file quite a bit as it's easily parsed.

Troubleshooting & further modification
--------------------------------------
If the scripts demand for files that I didn't explain how to obtain, let me know. I might have missed a file or two, it's--admittedly--a pretty convoluted pipeline.

The scripts ``create_go_annots_sprot_trembl.py`` and ``create_top_hit_sprot_trembl_nr.py`` have quite a few optional parameters to it. Check how they're used by calling the script with the ``-h`` flag, or read the source code.
