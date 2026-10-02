#!/usr/bin/env zsh

echo cd /Users/gmatz/Code/unipd/verlet3

mols=(0.5 0.75 1.0 1.25 1.5)
mols=(2)

for i in $(seq $#mols); do
mol_name="reactive-${mols[$i]}.xyz"
mol_name="gemini-${mols[$i]}.xyz"
mol_num=$((i - 1))
echo "#$mol_name"
echo mol delete $mol_num
echo "set molid [mol new {${mol_name}} type {xyz} first 0 last -1 step 1 waitfor -1]"

echo animate goto $mol_num
echo menu graphics on
echo mol modstyle 0 \$molid VDW 0.500000 12.000000
echo mol color Name
echo mol representation VDW 0.500000 12.000000
echo mol selection all
echo mol material Opaque
echo mol modrep 0 \$molid
echo menu graphics off
done
