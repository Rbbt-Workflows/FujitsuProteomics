require 'rbbt-util'
require 'rbbt/workflow'

Misc.add_libdir if __FILE__ == $0

#require 'rbbt/sources/FujitsuProteomics'

Workflow.require_workflow "Proteomics"

module FujitsuProteomics
  extend Workflow

  def self.organism
    "Hsa/feb2023"
  end

  def self.uni2ensp(uniprot)
    @@uni2ensp ||= Organism.protein_identifiers(organism).index :target => "Ensembl Protein ID", :persist => true
    @@uni2ensp2 ||= UniProt.identifiers[organism.split("/").first].index :target => "Ensembl Protein ID", :persist => true

    @@uni2ensp[uniprot] || @@uni2ensp2[uniprot]
  end

  input :variants_file, :file, "File with variants in extended format", nil, :required => true
  task :parse_mutated_isoform => :tsv do |variants_file|
    tsv = TSV.open variants_file, :type => :list, :key_field => 3, :namespace => FujitsuProteomics.organism
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
    tsv.column("Mutated Isoform").values.flatten.compact
  end


  dep :mutated_isoforms
  dep_task :annotate, Proteomics, :mi_wizard, :mutated_isoforms => :mutated_isoforms, :organism => FujitsuProteomics.organism

  dep :annotate, :compute => :produce
  task :annotate_file => :tsv do
    tsv = step(:parse_mutated_isoform).path.tsv :type => :double
    annotations = step(:annotate).load
    tsv.attach annotations
  end

  #{{{ DNA

  input :dna_variants_file, :file, "Clinical variants in DNA coordinates", nil, :required => true
  task :parse_genomic_variants => :array do |variants|
    TSV.traverse variants, :type => :array, :into => :stream do |line|
      chr, pos, ref, alt = line.split("\t")
      pos, muts = Misc.correct_vcf_mutation(pos, ref, alt)
      vars = muts.collect{|m| [chr, pos, m] * ":" }
      vars.extend MultipleResult
      vars
    end
  end

  dep :parse_genomic_variants
  dep_task :annotate_dna, Proteomics, :wizard, :mutations => :parse_genomic_variants, :organism => FujitsuProteomics.organism


end

#require 'FujitsuProteomics/tasks/basic.rb'

#require 'rbbt/knowledge_base/FujitsuProteomics'
#require 'rbbt/entity/FujitsuProteomics'

