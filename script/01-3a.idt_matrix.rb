
Fa1, Boudir, Fmat, min_aln_len, min_idt  = ARGV
Min_aln_len = min_aln_len.to_i
Min_idt     = min_idt.to_f

## parse labels and the order
labs = []
IO.read(Fa1).split(/^>/)[1..-1].each{ |ent|
  lab, *_seq = ent.split("\n")
  labs << lab.split(/\s+/)[0]
}

## parse blast output
mat = Hash.new{ |h, i| h[i] = {} }
Dir["#{Boudir}/*.blastp.out"].sort_by{ |fin| File.basename(fin).split(".")[0..-2]*"." }.each{ |fin|
  IO.readlines(fin).each{ |l|
    a = l.chomp.split("\t")
    que, sub, idt, aln_len = a.values_at(0, 1, 2, 3)

    next if aln_len.to_i < Min_aln_len
    next if idt.to_f < Min_idt

    mat[que][sub] = idt.to_f
  }
}

## make symmetric matrix
(0...labs.size).each{ |i|
  a = labs[i]
  (i...labs.size).each{ |j|
    b = labs[j]
    if mat[a][b] and mat[b][a]
      v = "%.6g" % ((mat[a][b] + mat[b][a]) / 2)
    elsif mat[a][b] 
      v = "%.6g" % mat[a][b] 
    elsif mat[b][a]
      v = "%.6g" % mat[b][a]
    else
      v = "0"
    end
    mat[a][b] = v
    mat[b][a] = v
  }
}

open(Fmat, "w"){ |fw|
  fw.puts ["", labs]*"\t"
  labs.each{ |lab1|
    vals = labs.map{ |lab2| mat[lab1][lab2] }
    fw.puts [lab1, vals]*"\t"
  }
}
