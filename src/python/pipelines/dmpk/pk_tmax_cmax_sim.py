#!/usr/bin/env python

# Copyright 2017 Informatics Matters Ltd.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
# http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.


# Estimation of absorption from t1/2 and Tmax after po administration (one-compartment)
# Based on original work by Amit Kumar Garg <a.garg@sygnaturediscovery.com>


import argparse, collections
import math   #math.log is the natural logarithm
from math import log10, floor
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from pipelines_utils import utils



### functions #########################################################

# TODO - move this function to utils
def round_sig(x, sig):
    return round(x, sig-int(floor(log10(abs(x))))-1)


def generatePlot(t_hf, t_hf_a, D, AUC, tn, quiet=False,
                 plot_height=4, plot_width=10, font_size=12, filename='cmax.png'):


    kel= math.log(2)/t_hf
    ka= math.log(2)/ t_hf_a
    Tmax=(math.log(ka)-math.log(kel))/(ka-kel)
    Cmax=math.exp(-kel*Tmax)*kel*AUC
    V_F=D/kel/AUC

    outputs = collections.OrderedDict()
    outputs['Tmax(hr)'] = round_sig(Tmax,3)
    outputs['Cmax(mg/L)'] = round_sig(Cmax,3)
    outputs['kel(hr-1)'] = round_sig(kel,3)
    outputs['ka(hr-1)'] = round_sig(ka,3)
    outputs['V/F(L)'] = round_sig(V_F,3)

    if not quiet:
        utils.log('------------------------------------------------------------------------------------------')
        utils.log('kel \t',kel)
        utils.log('ka \t',ka)
        utils.log('Tmax \t',Tmax)
        utils.log('Cmax \t',Cmax)
        utils.log('V_F \t',V_F)
        utils.log('------------------------------------------------------------------------------------------')

    b_time=[]
    c_cp=[]
    d_perc=[]
    for i in range(0,101):
        if(i==0):
            b_time.append(0)
        else:
            b_time.append(b_time[i-1]+tn/100)

        c_cp.append(ka*D/V_F/(ka-kel)*(math.exp(-kel*b_time[i])-math.exp(-ka*b_time[i])))

        d_perc.append(100-100*math.exp(-ka*b_time[i]))

    #print(b_time[100],c_cp[100],d_perc[100])

    # Creating the visualisation
    plt.figure(figsize=(plot_width,plot_height))
    plt.subplot(1, 2, 1)

    plt.plot(b_time,c_cp,linewidth=2,linestyle='dashed',color='coral')  #Plotting the observed data
    plt.xlabel('Time (h)',fontsize=font_size)
    plt.ylabel('Cp(mg/L',fontsize=font_size)
    plt.title('cp Vs Time',color='coral',fontsize=font_size)
    plt.grid(True)
    #plt.yscale('log')   #Change the Y sclae to logscale

    plt.subplot(1, 2, 2)
    plt.plot(b_time,d_perc,linewidth=2,linestyle='dashed')  #Plotting the observed data
    plt.xlabel('Time (h)',fontsize=font_size)
    plt.ylabel('% Absorbed',fontsize=font_size)
    plt.title('%Absorbed Vs Time',color='dodgerblue',fontsize=font_size)
    plt.grid(True)

    # Fine-tune figure; make subplots farther from each other.
    # refine layout to better support different sizes

    plt.savefig(filename)

    return outputs

### start main execution ##############################################

def main():

    ### command line args definitions ##################################

    parser = argparse.ArgumentParser(description='Tmax/Cmax simulation')
    parser.add_argument('--half-life', type=float, required=True, help='half life (hours)')
    parser.add_argument('--absorption', type=float, required=True, help='half life absorption (hours)')
    parser.add_argument('--dose', type=float, required=True, help='initial dose (mg)')
    parser.add_argument('--auc', type=float, required=True, help='AUC (mg/L*hr)')
    parser.add_argument('--time', type=float, required=True, help='time (h)')

    parser.add_argument('--plot-height', type=int, default=4, help='plot height')
    parser.add_argument('--plot-width', type=int, default=10, help='plot width')
    parser.add_argument('--font-size', type=int, default=12, help='font size')

    parser.add_argument('-o', '--output', type=str, default='output', help='output file base name')

    parser.add_argument('-q', '--quiet', action='store_true', help='Quiet mode')

    args = parser.parse_args()
    utils.log("Tmax/Cmax simulation Args: ", args)

    ### execute #######################################################

    outputs = generatePlot(args.half_life, args.absorption, args.dose, args.auc, args.time,
                 plot_width=args.plot_width, plot_height=args.plot_height, font_size=args.font_size,
                 filename=args.output + '.png')

    status_str = ", ".join(map(lambda e: e + ": " + str(outputs[e]), outputs))
    utils.write_metrics(args.output, {'__StatusMessage__':status_str, 'DMPK.Syg.TmaxCmax':1})


if __name__ == "__main__":
    main()