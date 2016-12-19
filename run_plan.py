import pickle as pkl

with open('plan.pickle','rb') as inpf:
  plan=pkl.load(inpf)

plan.nextstep()
