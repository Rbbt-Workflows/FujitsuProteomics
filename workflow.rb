require 'rbbt-util'
require 'rbbt/workflow'

Misc.add_libdir if __FILE__ == $0

#require 'rbbt/sources/FujitsuProteomics'

Workflow.require_workflow "Proteomics"

module FujitsuProteomics
  extend Workflow

  def self.organism
    "Hsa/may2017"
  end

  def self.uni2ensp(uniprot)
    @@uni2ensp ||= Organism.protein_identifiers(organism).index :target => "Ensembl Protein ID", :persist => true
    @@uni2ensp2 ||= UniProt.identifiers[organism.split("/").first].index :target => "Ensembl Protein ID", :persist => true

    @@uni2ensp[uniprot] || @@uni2ensp2[uniprot]
  end

  input :variants_file, :file, "File with variants in extended format", nil, :required => true
  task :parse_mutated_isoform => :tsv do |variants_file|
    tsv = TSV.open variants_file, :type => :list, :key_field => 3, :organism => FujitsuProteomics.organism
    tsv.key_field = "Variant"
    tsv.fields = %w(Gene HGVS AA Chromosome Start End UniProt)
    tsv.add_field "Mutated Isoform" do |key,values|
      gene, hgvs, aa, chr, start, eend, uniprot = values
      ensp = FujitsuProteomics.uni2ensp(uniprot)
      change = Misc.translate_prot_mutation_hgvs2rbbt(aa.dup)
      [ensp, change] * ":"
    end

    tsv
  end

  dep :parse_mutated_isoform
  task :mutated_isoforms => :array do
    tsv = step(:parse_mutated_isoform).load
    tsv.column("Mutated Isoform").values.flatten
  end


  dep :mutated_isoforms
  dep_task :annotate, Proteomics, :mi_wizard, :mutated_isoforms => :mutated_isoforms, :organism => FujitsuProteomics.organism

  dep :annotate
  task :annotate_file => :tsv do
    tsv = step(:parse_mutated_isoform).path.tsv :type => :double
    annotations = step(:annotate).load
    tsv.attach annotations
  end


end

#require 'FujitsuProteomics/tasks/basic.rb'

#require 'rbbt/knowledge_base/FujitsuProteomics'
#require 'rbbt/entity/FujitsuProteomics'

