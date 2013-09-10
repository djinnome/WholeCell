import gurobipy
def print_rxns( newRxns, oldRxns ):
    ratios = {}
    print "Rxn\tNewFlux\tOldFlux\tMultiplier"
    for rxn in sorted(newRxns):
        print "%s\t%f\t%f\t%f" % (rxn, newRxns[rxn], oldRxns[rxn], newRxns[rxn]/oldRxns[rxn])
        ratios[rxn] = newRxns[rxn]/oldRxns[rxn]
    return ratios

rxns = ['Tmk', 'Fba', 'MetK', 'Pgi', 'AceE', 'Pta', 'CmkA2', 'AckA', 'Eno', 'TpiA', 'Growth']
lpfile = '/msc/neurospora/src/WholeCell/data/metabolism.full.mmol-gDCW-h.gurobi.lp'
m = gurobipy.read(lpfile)
m.optimize()

# Store WT flux predictions into a dictionary
wtFlux = dict([(rxn.VarName, rxn.X) for rxn in m.getVars()])

# Store mutant fluxes predictions into a dictionary
ffp = ffparser.FFParser(
mutantFlux = dict([(rxn['Rxn'], float(rxn['Flux'])) for rxn in ffp.parse('../data/gold_fluxes.txt')[1:]])

#  f(x) = Lx + x^T*Q*x
# minimize the perturbed reaction fluxes subject to the constraint that growth rate is 3/4 wild type.
moma = gurobipy.quicksum([-mutantFlux[rxn]*m.getVarByName(rxn) for rxn in mutantFlux]) + gurobipy.quicksum([m.getVarByName(rxn)*m.getVarByName(rxn) for rxn in mutantFlux])




m.setObjective( moma , gurobipy.GRB.MINIMIZE )
m.update()
m.optimize()

# Store perturbed flux predictions into a dictionary
newRxns = dict([(rxn, m.getVarByName(rxn).X) for rxn in rxns])

# Print them out
ratios = print_rxns( newRxns, oldRxns )
print ratios
del m

# Confirm that setting the flux constraints to newRxns really does slow growth to 3/4  WT
m = gurobipy.read(lpfile)
for rxn in rxns:
    m.getVarByName(rxn).ub = newRxns[rxn]
m.update()
m.optimize()
