
require 'parallel'

Fa1, Boudir, Fmat, min_aln_len, min_idt, ncpus  = ARGV
Min_aln_len = min_aln_len.to_i
Min_idt     = min_idt.to_f
Ncpus       = ncpus.to_i

## parse labels and the order
labs = []
lab2idx = {}
IO.read(Fa1).split(/^>/)[1..-1].each.with_index{ |ent, idx|
  lab, *_seq = ent.split("\n")
  lab = lab.split(/\s+/)[0]

  labs << lab
  lab2idx[lab] = idx
}

## parse blast output
fins = Dir["#{Boudir}/*.blastp.out"].sort_by{ |fin| File.basename(fin).split(".")[0..-2]*"." }

### [!!!] blastp output should be split by query
mat = Parallel.map(fins, in_processes: Ncpus, progress: "parse blastp outputs"){ |fin|
  hit = Hash.new{ |h, i| h[i] = [] }
  _que = ""
  out = Array.new(labs.size) ### score array for each query, sized number of queries

  IO.readlines(fin).each{ |l|
    a = l.chomp.split("\t")
    que, sub, idt, aln_len  = a.values_at(0, 1, 2, 3)

    raise("blastp output should be split by query") if que != _que and _que != ""
    _que = que

    next if aln_len.to_i < Min_aln_len
    next if idt.to_f < Min_idt

    hit[sub] << a
  }

  hit.each{ |sub, as|
    pos2score = Hash.new(0)

    as.each{ |a|
      # qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore
      score, qst, qen    = a.values_at(11, 6, 7).map(&:to_i)
      qlen = qen - qst + 1

      new_scr = score.to_f / qlen.to_i ### score per query position

      (qst..qen).each{ |i|
        pos2score[i] = [pos2score[i], new_scr].max
      }
    }

    sum = 0
    pos2score.each{ |pos, scr|
      sum += scr
    }

    idx = lab2idx[sub]
    out[idx] = sum
  }


  ### store order of query
  _idx = lab2idx[_que]

  [out, _idx]
}

### sort out by query
scr = []
mat.each{ |out, _idx|
  scr[_idx] = out
}
p scr.size
(0...labs.size).each{ |i|
  begin
  p scr[i].size
  rescue
    p "error at #{i}"
    p scr[i]
    p labs[i]
  raise
  end
}

## make symmetric matrix
sim = []
labs.size.times do sim << Array.new(labs.size) end

(0...labs.size).each{ |i|
  (i...labs.size).each{ |j|
    a = scr[i][j]
    b = scr[j][i]
    _a = scr[i][i]
    _b = scr[j][j]

    if i == j
      v = "1"
    elsif a and b
      scr_a = a / _a
      scr_b = b / _b
      avg = (scr_a + scr_b) / 2
      v = "%.6g" % avg
    elsif a
      scr_a = a / _a
      v = "%.6g" % scr_a
    elsif b
      scr_b = b / _b
      v = "%.6g" % scr_b
    else
      v = "0"
    end
    sim[i][j] = v
    sim[j][i] = v
  }
}

## make output
open(Fmat, "w"){ |fw|
  fw.puts ["", labs]*"\t"
  (0...labs.size).each{ |i|
    lab1 = labs[i]
    fw.puts [lab1, sim[i]]*"\t"
  }
}
