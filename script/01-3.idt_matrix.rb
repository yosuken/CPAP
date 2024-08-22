
Fa1, Bout, Fmat, Min_aln_len = ARGV

## parse labels
labs = []
IO.read(Fa1).split(/^>/)[1..-1].each{ |ent|
  lab, *seq = ent.split("\n")
  labs << lab.split(/\s+/)[0]
}

## parse blast output
mat = Hash.new{ |h, i| h[i] = {} }
IO.readlines(Bout).each{ |l|
  a = l.chomp.split("\t")
  que, sub, idt, aln_len = a.values_at(0, 1, 2, 3)

  next if aln_len.to_i < Min_aln_len.to_i

  mat[que][sub] = idt.to_f
}

## make symmetric matrix
(0...labs.size).each{ |i|
  a = labs[i]
  (i...labs.size).each{ |j|
    b = labs[j]
    if mat[a][b] and mat[b][a]
      v = "%.2f" % ((mat[a][b] + mat[b][a]) / 2)
    elsif mat[a][b] 
      v = "%.2f" % mat[a][b] 
    elsif mat[b][a]
      v = "%.2f" % mat[b][a]
    else
      v = "0.00"
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
