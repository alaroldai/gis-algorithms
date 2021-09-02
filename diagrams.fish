# Assumes the data files are in gen/

for f in gen/*.dot
  dot -Tpng $f > $f.png
end

for f in gen/*.plt
  gnuplot -e "set terminal png size 800,800; set output '$f.png'" $f
end