
library(reticulate)
library(glue)

py_run_string(glue(
  "
import torch
import torch.nn as nn
import numpy as np
torch.set_default_tensor_type('torch.FloatTensor')
class GMMNet(nn.Module):
    def __init__(self,n_feature):
        super(GMMNet, self).__init__()
        self.fc = nn.Linear(n_feature,1)
        self.sig = nn.Sigmoid()
    def forward(self,UKBB_pop,theta_UKBB_GPC,study_info,colname_UKBB,var_SNP,var_GPC,C_):
        #UKBB_pop = torch.FloatTensor(UKBB_pop)
        #beta = torch.FloatTensor(beta)
        if isinstance(theta_UKBB_GPC,float):
          e10 =  torch.FloatTensor([theta_UKBB_GPC])
        else:
          e10 =  torch.FloatTensor([theta_UKBB_GPC[0]])
        e1 = self.fc(UKBB_pop[:,2:UKBB_pop.shape[1]])+e10
        e1 = self.sig(e1)
        e1 = e1.squeeze()
        u1 = torch.matmul((UKBB_pop[:,0] - e1),UKBB_pop[:,1:UKBB_pop.shape[1]])
        index_varSNP = [colname_UKBB.index(var_SNP[i]) for i in range(len(var_SNP))]
        u2_part1 = torch.matmul(e1,UKBB_pop[:,index_varSNP])
        if var_GPC is None:
          index_varGPC = []
        else:
          index_varGPC = [colname_UKBB.index(var_GPC[i]) for i in range(len(var_GPC))]
        
        u2_part2 = torch.empty(len(var_SNP))
        for snp_id in range(1,(len(var_SNP)+1)  ):#range(100,101):#

          index_snp_id = [0]+index_varGPC+[colname_UKBB.index('SNP'+str(snp_id))]
          u2_id = torch.matmul(UKBB_pop[:,index_snp_id],torch.FloatTensor(np.append(theta_UKBB_GPC,study_info[snp_id-1]['Coeff'])))
          u2_id = self.sig(u2_id)
          u2_part2[snp_id-1] =torch.dot(u2_id,UKBB_pop[:,colname_UKBB.index('SNP'+str(snp_id))])
        Un = torch.cat([u1,u2_part1-u2_part2],dim=0)
        Un = Un.squeeze()
        hatQ = torch.matmul(Un,torch.matmul(C_,Un))
        return hatQ/UKBB_pop.shape[0]**2



def torchoptimLBFGS(UKBB_pop,theta_UKBB_GPC,study_info,colname_UKBB,var_SNP,var_GPC,C_,initial_val,D,lam=1,lr=1,EPOCH=10,verbose=False):
    UKBB_pop = torch.FloatTensor(UKBB_pop)
    net = GMMNet(UKBB_pop.shape[1]-2)
    initial_val_tensor = torch.tensor(np.array(initial_val),dtype=torch.float)
    initial_val_tensor = initial_val_tensor[1:initial_val_tensor.shape[0]]
    initial_val_tensor = torch.nn.Parameter(initial_val_tensor.reshape(1,(UKBB_pop.shape[1]-2) ))
    net.fc.weight = initial_val_tensor
    net.fc.bias = torch.nn.Parameter(torch.tensor([0.]))
    net.fc.bias.requires_grad_(False)
    #for k in range(1,no_of_studies+1):
    #    study_info1[str(k)] = torch.FloatTensor(np.array(study_info1[str(k)]))
    C_ = torch.FloatTensor(np.array(C_))
    
    def null_loss(x):
        return x

    loss_metric = null_loss
    ## Initialize Optimizer and Learning Rate Scheduler
    learning_rate = lr
    optimizer = torch.optim.LBFGS(net.parameters(),max_iter=30)    
    ## scheduler = torch.optim.lr_scheduler.StepLR(optimizer, step_size=100, gamma=0.1)

    loss_collect = [100]
    count_early = 0
    eps_early = 10**(-10)


    def closure():
        outputs = net(UKBB_pop,theta_UKBB_GPC,study_info,colname_UKBB,var_SNP,var_GPC,C_)
        penalty_val = torch.tensor(0.)
        if lam != 0:
            pen_loc = np.where(np.diag(D)[1:D.shape[0]]==1)[0]
            
            penalty_val += (torch.norm(net.fc.weight[0,pen_loc],p=1))*lam
        optimizer.zero_grad() 
        loss = loss_metric(outputs)       
        loss += penalty_val
        loss.backward()
        return loss

    losses = optimizer.step(closure)
    i=0  
    #for i in range(N_LBFGS_STEPS_VALIDATION):
          
        #optimizer.zero_grad()
        #loss.backward()
        #optimizer.step()
        #scheduler.step()
    if (i+1)%100 == 0 or i == 0:
        if verbose:
            print('Epoch '+str(i+1)+' Loss: '+str(losses.item()))
    if verbose:
        print('Epoch '+str(i+1)+' Loss: '+str(losses.item()))
    loss_collect.append(losses.item())
    ## early stopping criteria

    #if np.abs(loss_collect[-1] - loss_collect[-2]) < eps_early:
    #    count_early += 1
    #if count_early >= 100:
    #    break
    
    if verbose:
        outputs = net(UKBB_pop,theta_UKBB_GPC,study_info,colname_UKBB,var_SNP,var_GPC,C_)
        penalty_val = torch.tensor(0.)
        if lam != 0:
            pen_loc = np.where(np.diag(D)[1:D.shape[0]]==1)[0]
            penalty_val += (torch.norm(net.fc.weight[0,pen_loc],p=1))*lam
        print('Final Loss: '+str((outputs+penalty_val).item()), flush = True)
    if isinstance(theta_UKBB_GPC,float):
      e10 =  [theta_UKBB_GPC]
    else:
      e10 =  [theta_UKBB_GPC[0]]
    return e10+[i.item() for i in net.fc.weight.squeeze()],np.array(loss_collect[1:len(loss_collect)])
")  
  ,local = FALSE, convert = TRUE)



# copy objects from the main python module into the specified R environment
py_main <- import_main(convert = TRUE)
py_main_dict <- py_get_attr(py_main, "__dict__")
py_name<- "torchoptimLBFGS"
Encoding(py_name) <- "UTF-8"

py_value <- py_main_dict[[py_name]]
py_envir<-globalenv()
assign(py_name, py_value, envir = py_envir) 

