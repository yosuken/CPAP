###
### CPAP (Clustering and Phylogenetic Analyzer of Proteins)
###
### Copyright: 2019-2025 (C) Yosuke Nishimura (nishimuray@jamstec.go.jp ; ynishimura@aori.u-tokyo.ac.jp)
### License: MIT license
###


# {{{ procedures
WriteBatch = lambda do |t, jobdir, outs|
  jdir = "#{jobdir}/#{t.name.split(":")[-1]}"; mkdir_p jdir unless File.directory?(jdir)
  jnum = outs.size

  if jnum > 0
    outs.each_slice(jnum).with_index(1){ |ls, idx| ## always 1 file
      open("#{jdir}/#{t.name.split(".")[1..-1]*"."}.sh", "w"){ |fjob|
        fjob.puts ls
      }
    }
  else
    open("#{jdir}/#{t.name.split(".")[1..-1]*"."}.sh", "w"){ |fjob| } ## clear job file
  end
end

RunBatch = lambda do |t, jobdir, ncpu, logdir|
  jdir = "#{jobdir}/#{t.name.split(":")[-1]}"
  ldir = "#{logdir}/#{t.name.split(":")[-1]}"; mkdir_p ldir unless File.directory?(ldir)

  Dir["#{jdir}/*.sh"].sort_by{ |fin| fin.split(".")[-1].to_i }.each{ |fin| ## always 1 or 0 file
    next if File.zero?(fin)
    sh "parallel --jobs #{ncpu} --joblog #{ldir}/parallel.log <#{fin}"
  }

  open("#{ldir}/exit", "w"){ |fw| fw.puts "exit at #{Time.now.strftime("%Y-%m-%d_%H:%M:%S")}" }
end

PrintStatus = lambda do |current, total, status, t|
  puts ""
  puts "\e[1;32m===== #{Time.now}\e[0m"
  puts "\e[1;32m===== step #{current} / #{total} (#{t.name}) -- #{status}\e[0m"
  puts ""
  $stdout.flush
end

CheckVersion = lambda do |commands|
  commands.each{ |command|
    str = case command
    when "ruby"
      %|ruby --version 2>&1|
    when "makeblastdb"
      %|makeblastdb -version 2>&1|
    when "blastp"
      %|blastp -version 2>&1|
    when "R"
      %|LANG=C R --version 2>&1|
    when "gplots"
      %|LANG=C R --quiet --no-save --no-restore -e "packageVersion('gplots')" 2>&1|
    when "phylogram"
      %|LANG=C R --quiet --no-save --no-restore -e "packageVersion('phylogram')" 2>&1|
    when "dendextend"
      %|LANG=C R --quiet --no-save --no-restore -e "packageVersion('dendextend')" 2>&1|
    when "ape"
      %|LANG=C R --quiet --no-save --no-restore -e "packageVersion('ape')" 2>&1|
    when "phangorn"
      %|LANG=C R --quiet --no-save --no-restore -e "packageVersion('phangorn')" 2>&1|
    end
    puts ""
    puts "\e[1;32m===== check version: #{command}\e[0m"
    puts ""
    puts "$ #{str}"
    ### run
    puts `#{str}`
    ### flush
    $stdout.flush
  }
end
# }}} procedures


# {{{ task controller
task :default do
  ## argument file path
  Odir       = ENV["dir"]            ## base output directory
  Fin        = ENV["fin"]            ## input protein fasta file

  ### check 'parallel' avilability
  require 'parallel'

  ## phylogenetic tree file (newick or nexus)
  if ENV["fphy"] != ""
    Fphy  = ENV["fphy"]
  elsif ENV["fphy_as_chronogram"] != ""
    Fphy  = ENV["fphy_as_chronogram"]
  else
    Fphy  = ""
  end
  Use_chrono = ENV["fphy_as_chronogram"] != "" ? true : false  ## flag if to use chronogram

  ## file path to be generated
  Jobdir     = "#{Odir}/batch"
  Logdir     = "#{Odir}/log/tasks"
  Bdbdir     = "#{Odir}/blastdb"      ## dir for blastdb
  Boudir     = "#{Odir}/blastp"       ## dir for blast output
  Resdir     = "#{Odir}/result"       ## dir for result
  Spldir     = "#{Odir}/input_split"  ## dir for input_split
  Fa1        = "#{Odir}/input.faa"    ## copy of Fin
  Fa2        = "#{Bdbdir}/input.faa"  ## copy of Fin (as symlink of Fa1)

  ## measure
  Measure    = ENV["measure"]               
  case Measure
  when "identity"  ## %identity of (default)
    lab      = "idt"
    Fmat     = "#{Resdir}/#{lab}.tsv"     ## identity matrix
    Fpdf     = "#{Resdir}/#{lab}-heat"    ## heatmap pdf
    Fdnd     = "#{Resdir}/dendrogram"     ## dendrogram
    FmatR    = "#{Resdir}/#{lab}.ordered" ## identity matrix in the order of the dendrogram

    ## tasks
    tasks    = %w|01-1.makeblastdb 01-2.blastp 01-3a.idt_matrix 01-4.heatmap 01-5.reordered_idt_matrix|
  when "sim-score" ## use Sg like similarity score
    ### [TODO] implement calculation of sim-score
    lab      = "sim"
    Fmat     = "#{Resdir}/#{lab}.tsv"     ## identity matrix
    Fpdf     = "#{Resdir}/#{lab}-heat"    ## heatmap pdf
    Fdnd     = "#{Resdir}/dendrogram"     ## dendrogram
    FmatR    = "#{Resdir}/#{lab}.ordered" ## identity matrix in the order of the dendrogram

    ## tasks
    tasks    = %w|01-1.makeblastdb 01-2.blastp 01-3b.sim_matrix 01-4.heatmap 01-5.reordered_idt_matrix|
    # tasks    = %w|01-3b.sim_matrix|
  else ## not defined
    raise("`--measure #{Measure}': does not defined")
  end


  ## params
  # computing
  Ncpus           = ENV["ncpus"]        ## default: 4
  # BLASTp
  DBsize          = ENV["dbsize"]       ## default: 100,000,000
  Matrix          = ENV["matrix"]       ## default: BLOSUM62
  Evalue          = ENV["evalue"]       ## default: 0.01
  Min_aln_len     = ENV["min_aln_len"]  ## default: 40
  Min_idt         = ENV["min_idt"]      ## default: 20
  Max_target_seqs = 1_000_000
  # hclust clustering
  Clust_method    = ENV["clust_method"] ## default: average (option: ward.D2, ...)
  Clust_methods   = %w|average ward.D ward.D2 single complete mcquitty median centroid|

  unless (Clust_methods + %w|ALL|).include?(Clust_method)
    raise("--clust-method #{Clust_method}: unknown clustering method. aborting...")
  end

  ### check version
  commands = %w|blastp makeblastdb R dendextend gplots ruby|
  CheckVersion.call(commands)

  ### run
  NumStep = tasks.size
  tasks.each.with_index(1){ |task, idx|
    Rake::Task[task].invoke(idx)
  }
end
# }}}


# {{{ "01-1.makeblastdb"
desc "01-1.makeblastdb"
task "01-1.makeblastdb", ["step"] do |t, args|
  PrintStatus.call(args.step, NumStep, "START", t)
  mkdir_p Odir
  mkdir_p Bdbdir
  log  = "#{Bdbdir}/makeblastdb.log"

  sh "cp #{Fin} #{Fa1}" unless File.exist?(Fa1)  ## copy input fasta to output dir
  sh "pushd #{Bdbdir} ; ln -s ../#{File.basename(Fa1)} ; popd" unless File.exist?(Fa2)
  sh "pushd #{Bdbdir} ; makeblastdb -dbtype prot -in #{File.basename(Fa2)} -out #{File.basename(Fa2)} 2>#{File.basename(log)} ; popd" unless File.exist?(log)

  ### split input fasta (one file for one seq)
  sh "mkdir -p #{Spldir}" unless File.directory?(Spldir)
  IO.read(Fa1).split(/^>/)[1..-1].each{ |ent|
    lab, *_seq = ent.split("\n")
    gid = lab.split(/\s+/)[0]
    fou = "#{Spldir}/#{gid}.faa"
    open(fou, "w"){ |fw| fw.puts ">#{ent}" } unless File.exist?(fou)
  }
end
# }}}


# {{{ "01-2.blastp"
desc "01-2.blastp"
task "01-2.blastp", ["step"] do |t, args|
  PrintStatus.call(args.step, NumStep, "START", t)
  outs   = []

  mkdir_p Boudir
  outfmt = '-outfmt "6 std qlen slen"' ### equal to '-outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen"'

  Dir["#{Spldir}/*.faa"].sort_by{ |fin| File.basename(fin).split(".")[0..-2]*"." }.each{ |fin|
    gid  = File.basename(fin).split(".")[0..-2]*"."
    bout = "#{Boudir}/#{gid}.blastp.out"
    blog = "#{Boudir}/#{gid}.blastp.log"

    outs << "blastp -num_threads 1 -matrix #{Matrix} -evalue #{Evalue} -dbsize #{DBsize} -max_target_seqs #{Max_target_seqs} -db #{Fa2} -query #{fin} -out #{bout} #{outfmt} 2>#{blog}" unless File.exist?(blog)
  }

  WriteBatch.call(t, Jobdir, outs)
  RunBatch.call(t, Jobdir, Ncpus, Logdir)
end
# }}}


# {{{ "01-3a.idt_matrix"
desc "01-3a.idt_matrix"
task "01-3a.idt_matrix", ["step"] do |t, args|
  PrintStatus.call(args.step, NumStep, "START", t)
  mkdir_p Resdir unless File.directory?(Resdir)
  script = "#{File.dirname(__FILE__)}/script/#{t.name}.rb"

  sh "ruby #{script} #{Fa1} #{Boudir} #{Fmat} #{Min_aln_len} #{Min_idt}" unless File.exist?(Fmat)
end
# }}}


# {{{ "01-3b.sim_matrix"
desc "01-3b.sim_matrix"
task "01-3b.sim_matrix", ["step"] do |t, args|
  PrintStatus.call(args.step, NumStep, "START", t)
  mkdir_p Resdir unless File.directory?(Resdir)
  script = "#{File.dirname(__FILE__)}/script/#{t.name}.rb"

  sh "ruby #{script} #{Fa1} #{Boudir} #{Fmat} #{Min_aln_len} #{Min_idt} #{Ncpus}" unless File.exist?(Fmat)
end
# }}}


# {{{ "01-4.heatmap"
desc "01-4.heatmap"
task "01-4.heatmap", ["step"] do |t, args|
  PrintStatus.call(args.step, NumStep, "START", t)
  script = "#{File.dirname(__FILE__)}/script/#{t.name}.R"

  case Clust_method
  when "ALL"
    Clust_methods.each{ |cm|
      fpdf = "#{Fpdf}.#{cm}.pdf"
      fdnd = "#{Fdnd}.#{cm}.newick"

      cmd  = "LANG=C Rscript #{script} #{Fmat} #{fpdf} #{fdnd} #{cm}"
      cmd  = "#{cmd} #{Fphy}"    if Fphy != "" ## when --fphy or --fphy_as_chronogram is given
      cmd  = "#{cmd} chronogram" if Use_chrono ## when --fphy_as_chronogram is given

      sh cmd unless File.exist?(Fdnd)
    }
  else
    cm   = Clust_method
    fpdf = "#{Fpdf}.#{cm}.pdf"
    fdnd = "#{Fdnd}.#{cm}.newick"

    cmd  = "LANG=C Rscript #{script} #{Fmat} #{fpdf} #{fdnd} #{cm}"
    cmd  = "#{cmd} #{Fphy}"    if Fphy != "" ## when --fphy or --fphy_as_chronogram is given
    cmd  = "#{cmd} chronogram" if Use_chrono ## when --fphy_as_chronogram is given

    sh cmd unless File.exist?(Fdnd)
  end
end
# }}}


# {{{ "01-5.reordered_idt_matrix"
desc "01-5.reordered_idt_matrix"
task "01-5.reordered_idt_matrix", ["step"] do |t, args|
  PrintStatus.call(args.step, NumStep, "START", t)
  mkdir_p Resdir
  script   = "#{File.dirname(__FILE__)}/script/#{t.name}.rb"

  case Clust_method
  when "ALL"
    Clust_methods.each{ |cm|
      fdnd  = "#{Fdnd}.#{cm}.newick"
      fmatR = "#{FmatR}.#{cm}.tsv"

      sh "ruby #{script} #{Fmat} #{fdnd} #{fmatR}" unless File.exist?(fmatR)
    }
  else
    cm    = Clust_method
    fdnd  = "#{Fdnd}.#{cm}.newick"
    fmatR = "#{FmatR}.#{cm}.tsv"

    sh "ruby #{script} #{Fmat} #{fdnd} #{fmatR}" unless File.exist?(fmatR)
  end
end
# }}}
