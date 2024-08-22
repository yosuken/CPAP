
Fmat, Fdnd, FmatR = ARGV

## parse dendrogram
# (('a':0,'b':0):190.2486,('c':164.5815,('d':130.3324,('e':117.1488,('G':107.2966,('f':99.37408,('g':94.174,('h':17.81668,'i':17.81668):76.35732):5.200086):7.922551):9.852154):13.18362):34.24904):25.66714);
gids = {}
idx  = 0
str = IO.readlines(Fdnd)[0]
str.split(/[\(\),;\s]+/).each{ |i|
  e = i.split(":")[0]
  if e and e != ""
    gid = e.gsub("'", "")

    gids[gid] = idx ; idx += 1
  end
}

## parse matrix
mat = []
n = gids.size
n.times do mat << Array.new(n) end

head, *body = IO.readlines(Fmat)

lab, *_gids = head.chomp.split("\t")
c = {} ### new_idx --> old_idx
_gids.each.with_index{ |gid, old|
  new = gids[gid] 
  c[new] = old
}

body.each.with_index{ |l, i|
  gid, *a = l.chomp.split("\t", -1)

  a.each.with_index{ |v, j|
    mat[i][j] = v
  }
}

open(FmatR, "w"){ |fw|
  fw.puts [lab, gids.keys]*"\t"

  (0..n-1).each{ |i|
    a = [gids.keys[i]]
    (0..n-1).each{ |j|
      a << mat[c[i]][c[j]]
    }

    fw.puts a*"\t"
  }
}
