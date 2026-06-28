trainlabel = [-1;1];
traindata = [2,1;3,4];

model = train(trainlabel,sparse(traindata));

int i = 1;