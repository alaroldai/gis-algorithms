# Assumes the data files are in gen/

for f in gen/*.dot
  set out img/(basename $f)
  dot -Tpng $f > $out.png
end

for f in gen/*.plt
  set out img/(basename $f)
  gnuplot -e "set terminal png size 800,800; set output '$out.png'" $f
end