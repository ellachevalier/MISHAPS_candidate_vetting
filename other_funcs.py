#!/usr/bin/env python
# coding: utf-8
# %%

# %%
from PyAstronomy.pyasl import foldAt
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
from matplotlib import gridspec
import matplotlib.image as mpimg
from pypdf import PdfMerger
from tabulate import tabulate
from texttable import Texttable
import latextable
import aspose.pdf as ap
from fpdf import FPDF
import fpdf
from PyPDF2 import PdfWriter, PdfReader
import json

def make_chi2_csv(initial_best_periods, best_periods, field, chip, star_id):
    filename1='MISHAPs_'+str(field)+'_'+str(chip)+'_'+str(star_id)+'_initial_periods_chi2.csv'
    filename2='MISHAPs_'+str(field)+'_'+str(chip)+'_'+str(star_id)+'_final_periods_chi2.csv'
    initial_best_periods.to_csv(os.path.join('period_fit_results',filename1), index=False)
    best_periods.to_csv(os.path.join('period_fit_results',filename2), index=False)
    
def make_json_file(chains, final_periods, field, chip, star_id):
    filename='MISHAPs_'+str(field)+'_'+str(chip)+'_'+str(star_id)+'_mcmc_chains.json'
    json_object = json.dumps(chains, indent=4)
 
    # Writing to sample.json
    with open(os.path.join('chains',filename), "w") as outfile:
        outfile.write(json_object)   

def make_final_pdf(contents, field, chip, star_id, p0_guess, even_guesses, odd_guesses, all_guesses, final_periods, densities_r, densities_z, best_periods, transit_nights, remove_variability=False):
    starid='MISHAPS_'+str(field)+'_'+str(chip)+'_'+str(star_id)
    with open("tables.tex", "w") as file:
        file.write('\documentclass[8pt]{article}\n')
        file.write('\\begin{document}\n')
        file.write('\\addtolength{\oddsidemargin}{-1.5in}\n')
        file.write('\\addtolength{\evensidemargin}{-1.5in}\n')
        for content in contents:
            file.write('\\vspace{-6in}\n')
            file.write(content)
        file.write('\n')
        file.write('\end{document}')
    

    options = ap.TeXLoadOptions()

    document = ap.Document("tables.tex" , options)
    # Convert Latex to PDF
    document.save("tables.pdf")
    
    
    p0_content='p0 guess: '+str(np.array(p0_guess, dtype=object))
    even_content='even guess: '+str(even_guesses)
    odd_content='odd guess: '+str(odd_guesses)
    all_content='all guess: '+str(all_guesses)
    pdf = FPDF()
    pdf.add_page()
    pdf.set_font('arial', 'B', 16.0)
    pdf.set_top_margin(8.0)
    pdf.ln(10.0)
    pdf.cell(ln=1, h=5.0, align='L', w=0, txt='Field: ' +str(field))
    pdf.cell(ln=1, h=5.0, align='L', w=0, txt='Chip: ' +str(chip))
    pdf.cell(ln=1, h=5.0, align='L', w=0, txt='Target: ' +str(star_id))
    pdf.set_font('arial', '', 6.0)
    pdf.ln(6.0)
    pdf.cell(ln=1, h=3.0, align='L', w=0, txt=p0_content, border=0)
    pdf.cell(ln=1, h=3.0, align='L', w=0, txt=even_content, border=0)
    pdf.cell(ln=1, h=3.0, align='L', w=0, txt=odd_content, border=0)
    pdf.cell(ln=1, h=3.0, align='L', w=0, txt=all_content, border=0)
    pdf.ln(6.0)
    pdf.set_font('arial', '', 12.0)
    pdf.cell(ln=10, h=3.0, align='L', w=0, txt='Densities from all folded data:', border=0)
    pdf.ln(2.0)
    pdf.set_font('arial', '', 10.0)
    
    for period in final_periods:
        pdf.cell(ln=1, h=3.0, align='L', w=0, txt='Period: '+str(period))
        pdf.ln(1.0) 
        pdf.cell(ln=1, h=3.0, align='L', w=0, txt='r-band: '+str(densities_r[period][0])+u" \u00B1 "+str(densities_r[period][1]), border=0)
        pdf.ln(1.0) 
        pdf.cell(ln=1, h=3.0, align='L', w=0, txt='z-band: '+str(densities_z[period][0])+u" \u00B1 "+str(densities_z[period][1]), border=0)
        pdf.ln(3.0)
        
    pdf.ln(6.0)
    pdf.cell(ln=1, h=3.0, align='L', w=0, txt='Best periods output:')
    pdf.ln(2.0)
    pdf.cell(ln=1, h=3.0, align='L', w=0, txt='First transit: '+str(transit_nights[0])+'    Last transit: '+str(transit_nights[-1]))
    pdf.ln(4.0)
    pdf.cell(ln=1, h=3.0, align='L', w=0, txt='Period'+ ' '*26+ 'x2_transit'+ ' '*14+'N')
    pdf.ln(2.0)
    
    if len(best_periods)<20:
        cutoff=len(best_periods)
    else:
        cutoff=20
    period_list=list(best_periods['Period'])
    x2_list=list(best_periods['x2_transit'])
    n_list=list(best_periods['N'])
    for idx in range(cutoff):
        string=str(round(period_list[idx],10))+' '*8+str(round(x2_list[idx],6))+' '*8+str(round(n_list[idx]))
        pdf.cell(ln=1, h=3.0, align='L', w=0, txt=string, border=0)
        pdf.ln(3.0)
        
    
    pdf.output('guesses.pdf', 'F')
    
    result_pdfs=[]
    for best_period in final_periods:
        result_pdfs.append(os.path.join('results', "results_"+str(starid)+"_period_"+str(round(best_period,4))+".pdf"))
    #pdfs = ['final_PDF.pdf']
    pdfs = ['guesses.pdf', 'tables.pdf']
    pdfs.append(os.path.join('figs','individual_transits_plot.pdf'))
    pdfs.extend(result_pdfs)
    remove_variability_pdfs = ['var_cut_lc_check.pdf', 'var_before_removal_check.pdf', 'variability_removal_check.pdf']
    if remove_variability==True:
        for file in remove_variability_pdfs:
            pdfs.append(os.path.join('remove_variability_results', file))
    period_pdfs = ['log_x2s_vs_log_periods_with_final_periods.pdf', 'deltax2_vs_period.pdf', 'total_plot_'+str(star_id)+'.pdf']
    for file in period_pdfs:
        pdfs.append(os.path.join('period_fit_results', file))
    
    fig_list=os.listdir('figs')
    sec_pic_list = [item for item in fig_list if 'secondary' in item]
    for file in sec_pic_list:
        pdfs.append(os.path.join('figs', file))
    #pdfs.append(os.path.join('figs', 'secondary_eclipse_fit_all_parameters.pdf'))
    #pdfs.append('total_plot_'+str(star_id)+'.pdf')
    secondary_pdfs=[]
    for best_period in final_periods:
        secondary_pdfs.append(os.path.join('results',"secondary_results"+"_period_"+str(round(best_period,4))+".pdf"))
    pdfs.extend(secondary_pdfs)
    merger = PdfMerger()
    for pdf in pdfs:
        merger.append(os.path.join(pdf))

    if os.path.exists('results')==False:
        os.mkdir('results')
    merger.write(os.path.join('results',"final_pdf_"+str(starid)+".pdf"))
    merger.close()

def make_latex_tables(results_batman_all, results_batman_even, results_batman_odd, results_mcmc_per_all, results_mcmc_per_even, results_mcmc_per_odd, results_mcmc_all, results_mcmc_even, results_mcmc_odd, best_period):
    all_data_list =[str(round(results_batman_all['t0'],8))+'$\pm$'+str(round(results_batman_all['t0_err'],8)),
              str(round(results_batman_all['rp_r'],8))+'$\pm$'+str(round(results_batman_all['rp_r_err'],8)),
              str(round(results_batman_all['rp_z'],8))+'$\pm$'+str(round(results_batman_all['rp_z_err'],8)),
              str(round(results_batman_all['a'],8))+'$\pm$'+str(round(results_batman_all['a_err'],8)),
              str(round(results_batman_all['inc'],8))+'$\pm$'+str(round(results_batman_all['inc_err'],8)),
              str(round(results_batman_all['C_r'],8))+'$\pm$'+str(round(results_batman_all['C_r_err'],8)),
              str(round(results_batman_all['C_z'],8))+'$\pm$'+str(round(results_batman_all['C_z_err'],8))]
    if results_batman_even is not None:
        even_data_list =[str(round(results_batman_even['t0'],8))+'$\pm$'+str(round(results_batman_even['t0_err'],8)),
                  str(round(results_batman_even['rp_r'],8))+'$\pm$'+str(round(results_batman_even['rp_r_err'],8)),
                  str(round(results_batman_even['rp_z'],8))+'$\pm$'+str(round(results_batman_even['rp_z_err'],8)),
                  str(round(results_batman_even['a'],8))+'$\pm$'+str(round(results_batman_even['a_err'],8)),
                  str(round(results_batman_even['inc'],8))+'$\pm$'+str(round(results_batman_even['inc_err'],8)),
                  str(round(results_batman_even['C_r'],8))+'$\pm$'+str(round(results_batman_even['C_r_err'],8)),
                  str(round(results_batman_even['C_z'],8))+'$\pm$'+str(round(results_batman_even['C_z_err'],8))]
    else:
        even_data_list = ['n/a', 'n/a', 'n/a', 'n/a', 'n/a', 'n/a', 'n/a']
    if results_batman_odd is not None:
        odd_data_list =[str(round(results_batman_odd['t0'],8))+'$\pm$'+str(round(results_batman_odd['t0_err'],8)),
                  str(round(results_batman_odd['rp_r'],8))+'$\pm$'+str(round(results_batman_odd['rp_r_err'],8)),
                  str(round(results_batman_odd['rp_z'],8))+'$\pm$'+str(round(results_batman_odd['rp_z_err'],8)),
                  str(round(results_batman_odd['a'],8))+'$\pm$'+str(round(results_batman_odd['a_err'],8)),
                  str(round(results_batman_odd['inc'],8))+'$\pm$'+str(round(results_batman_odd['inc_err'],8)),
                  str(round(results_batman_odd['C_r'],8))+'$\pm$'+str(round(results_batman_odd['C_r_err'],8)),
                  str(round(results_batman_odd['C_z'],8))+'$\pm$'+str(round(results_batman_odd['C_z_err'],8))]
    else:
        odd_data_list = ['n/a', 'n/a', 'n/a', 'n/a', 'n/a', 'n/a', 'n/a']
        
    df_batman = pd.DataFrame(dict(param=['t0', 'rpr', 'rpz', 'a', 'inc', 'Cr', 'Cz'],
        all_data=all_data_list, even_data=even_data_list, odd_data=odd_data_list))
    dict_new = {'param':'Param',
        'all_data': 'All',
        'even_data': 'Even',
        'odd_data': 'Odd'}
 
    # call rename () method
    df_batman.rename(columns=dict_new,
              inplace=True)     
    Row_list =[['Param', 'All', 'Even', 'Odd']]
    #Row_list =[['data', 't0', 'rpr', 'rpz', 'a', 'inc', 'Cr', 'Cz']]
             
    # Iterate over each row
    for index, rows in df_batman.iterrows():
            # Create list for the current row
        my_list =[rows.Param, rows.All, rows.Even, rows.Odd]

            # append the list to the final list
        Row_list.append(my_list)
    table1 = Texttable()
    table1.set_cols_align(["c"]*4)
    table1.set_deco(Texttable.HEADER | Texttable.VLINES)
    table1.add_rows(Row_list)
    content1=latextable.draw_latex(table1, caption='Batman parameters for P='+str(best_period))
    
    all_data_list_mcmc=[str(round(results_mcmc_per_all['50th'][0],8))+'$\pm$'+str(round(results_mcmc_all['t0_err'],8)),
                 str(round(results_mcmc_per_all['50th'][1],8))+'$\pm$'+str(round(results_mcmc_all['rp_r_err'],8)),
                 str(round(results_mcmc_per_all['50th'][2],8))+'$\pm$'+str(round(results_mcmc_all['rp_z_err'],8)),
                 str(round(results_mcmc_per_all['50th'][3],8))+'$\pm$'+str(round(results_mcmc_all['a_err'],8)),
                 str(round(results_mcmc_per_all['50th'][4],8))+'$\pm$'+str(round(results_mcmc_all['inc_err'],8)),
                 str(round(results_mcmc_per_all['50th'][5],8))+'$\pm$'+str(round(results_mcmc_all['C_r_err'],8)),
                 str(round(results_mcmc_per_all['50th'][6],8))+'$\pm$'+str(round(results_mcmc_all['C_z_err'],8))]
    if results_mcmc_even is not None:
        even_data_list_mcmc=[str(round(results_mcmc_per_even['50th'][0],8))+'$\pm$'+str(round(results_mcmc_even['t0_err'],8)),
                     str(round(results_mcmc_per_even['50th'][1],8))+'$\pm$'+str(round(results_mcmc_even['rp_r_err'],8)),
                     str(round(results_mcmc_per_even['50th'][2],8))+'$\pm$'+str(round(results_mcmc_even['rp_z_err'],8)),
                     str(round(results_mcmc_per_even['50th'][3],8))+'$\pm$'+str(round(results_mcmc_even['a_err'],8)),
                     str(round(results_mcmc_per_even['50th'][4],8))+'$\pm$'+str(round(results_mcmc_even['inc_err'],8)),
                     str(round(results_mcmc_per_even['50th'][5],8))+'$\pm$'+str(round(results_mcmc_even['C_r_err'],8)),
                     str(round(results_mcmc_per_even['50th'][6],8))+'$\pm$'+str(round(results_mcmc_even['C_z_err'],8))]
    else:
        even_data_list_mcmc=['n/a', 'n/a', 'n/a', 'n/a', 'n/a', 'n/a', 'n/a']
    if results_mcmc_odd is not None:
        odd_data_list_mcmc=[str(round(results_mcmc_per_odd['50th'][0],8))+'$\pm$'+str(round(results_mcmc_odd['t0_err'],8)),
                     str(round(results_mcmc_per_odd['50th'][1],8))+'$\pm$'+str(round(results_mcmc_odd['rp_r_err'],8)),
                     str(round(results_mcmc_per_odd['50th'][2],8))+'$\pm$'+str(round(results_mcmc_odd['rp_z_err'],8)),
                     str(round(results_mcmc_per_odd['50th'][3],8))+'$\pm$'+str(round(results_mcmc_odd['a_err'],8)),
                     str(round(results_mcmc_per_odd['50th'][4],8))+'$\pm$'+str(round(results_mcmc_odd['inc_err'],8)),
                     str(round(results_mcmc_per_odd['50th'][5],8))+'$\pm$'+str(round(results_mcmc_odd['C_r_err'],8)),
                     str(round(results_mcmc_per_odd['50th'][6],8))+'$\pm$'+str(round(results_mcmc_odd['C_z_err'],8))]
    else:
        odd_data_list_mcmc=['n/a', 'n/a', 'n/a', 'n/a', 'n/a', 'n/a', 'n/a']
    df_mcmc = pd.DataFrame(dict(param=['t0', 'rpr', 'rpz', 'a', 'inc', 'Cr', 'Cz'],
        all_data=all_data_list_mcmc, even_data=even_data_list_mcmc, odd_data=odd_data_list_mcmc))
    dict_new = {'param':'Param',
        'all_data': 'All',
        'even_data': 'Even',
        'odd_data': 'Odd'}
 
    # call rename () method
    df_mcmc.rename(columns=dict_new,
          inplace=True)    
    
    #Row_list =[['data', 't0', 'rpr', 'rpz', 'a', 'inc', 'Cr', 'Cz']]
    Row_list =[['Param', 'All', 'Even', 'Odd']]
    # Iterate over each row
    for index, rows in df_mcmc.iterrows():
            # Create list for the current row
        my_list =[rows.Param, rows.All, rows.Even, rows.Odd]

            # append the list to the final list
        Row_list.append(my_list)
    table2 = Texttable()
    table2.set_cols_align(["c"]*4)
    table2.set_deco(Texttable.HEADER | Texttable.VLINES)
    table2.add_rows(Row_list)
    content2=latextable.draw_latex(table2, caption='Mcmc parameters for P='+str(best_period))
    
    return [content1, content2]

def make_pdf(best_period, starid):
    fig_list=os.listdir('figs')
    #print(fig_list)
    even_pic_list = [item for item in fig_list if 'Even' in item and 'png' in item]
    odd_pic_list = [item for item in fig_list if 'Odd' in item and 'png' in item]
    all_pic_list = [item for item in fig_list if 'All' in item and 'png' in item]
    even_pic_list.sort()
    odd_pic_list.sort()
    all_pic_list.sort()
    #pic_list
    #fig_list
    #print(pic_list)
    
    model_all=mpimg.imread(os.path.join('figs', all_pic_list[0]))
    compare_all=mpimg.imread(os.path.join('figs', all_pic_list[1]))
    if len(even_pic_list)!=0:
        model_even=mpimg.imread(os.path.join('figs', even_pic_list[0]))
        compare_even=mpimg.imread(os.path.join('figs', even_pic_list[1]))
    if len(odd_pic_list)!=0:
        model_odd=mpimg.imread(os.path.join('figs', odd_pic_list[0]))
        compare_odd=mpimg.imread(os.path.join('figs', odd_pic_list[1]))
    
    fig=plt.figure(figsize=(14,8))
    fig.tight_layout()
    rows=2
    columns=3

    plt.subplot(rows,columns,1)
    plt.imshow(model_all)
    plt.axis('off')
    
    plt.subplot(rows,columns,4)
    plt.imshow(compare_all)
    plt.axis('off')

    if len(even_pic_list)!=0:
        plt.subplot(rows,columns,2)
        plt.imshow(model_even)
        plt.axis('off')

        plt.subplot(rows,columns,5)
        plt.imshow(compare_even)
        plt.axis('off')
        
    if len(odd_pic_list)!=0:
        plt.subplot(rows,columns,3)
        plt.imshow(model_odd)
        plt.axis('off')

        plt.subplot(rows,columns,6)
        plt.imshow(compare_odd)
        plt.axis('off')

    fig.subplots_adjust(wspace=0,hspace=0)

    fig.savefig(os.path.join('figs', 'models_and_rp_corners.pdf'))
    
    final_list=['comparison_of_rp.pdf', 'models_and_rp_corners.pdf']
    fig_list=os.listdir('figs')
    fig_list=[item for item in fig_list if ('walkers' in item or 'corner_plot' in item) and 'Individual' not in item]
    all_list=[item for item in fig_list if 'All' in item]
    even_list=[item for item in fig_list if 'Even' in item]
    odd_list=[item for item in fig_list if 'Odd' in item]
    final_list.extend(all_list)
    final_list.extend(even_list)
    final_list.extend(odd_list)
    final_list
    
    pdfs = final_list

    merger = PdfMerger()

    for pdf in pdfs:
        merger.append(os.path.join('figs', pdf))

    if os.path.exists('results')==False:
        os.mkdir('results')
    merger.write(os.path.join('results',"results_"+str(starid)+"_period_"+str(round(best_period,4))+".pdf"))
    merger.close()
    
    for item in even_pic_list:
        os.remove(os.path.join('figs', item))
    for item in odd_pic_list:
        os.remove(os.path.join('figs', item))
    for item in all_pic_list:
        os.remove(os.path.join('figs', item))
    #for f in os.listdir('figs'):
        #os.remove(os.path.join('figs', f))


def comparison_plot(results_mcmc_even, results_mcmc_odd, results_mcmc_all, results_batman_even, results_batman_odd, results_batman_all, period=None):
    
    if results_mcmc_even is not None:
        rp_r_even_batman=results_batman_even['rp_r']
        rp_r_even_err_batman=results_batman_even['rp_r_err']
        rp_z_even_batman=results_batman_even['rp_z']
        rp_z_even_err_batman=results_batman_even['rp_z_err']
        rp_r_even=results_mcmc_even['50th'][1]
        rp_r_err_even16=results_mcmc_even['16th'][1]
        rp_r_err_even84=results_mcmc_even['84th'][1]
        rp_z_even=results_mcmc_even['50th'][2]
        rp_z_err_even16=results_mcmc_even['16th'][2]
        rp_z_err_even84=results_mcmc_even['84th'][2]
        
        rp_r_even_med=results_mcmc_even['Median'][1]
        rp_z_even_med=results_mcmc_even['Median'][2]
        
    if results_mcmc_odd is not None:
        rp_r_odd_batman=results_batman_odd['rp_r']
        rp_r_odd_err_batman=results_batman_odd['rp_r_err']
        rp_z_odd_batman=results_batman_odd['rp_z']
        rp_z_odd_err_batman=results_batman_odd['rp_z_err']
        rp_r_odd=results_mcmc_odd['50th'][1]
        rp_r_err_odd16=results_mcmc_odd['16th'][1]
        rp_r_err_odd84=results_mcmc_odd['84th'][1]
        rp_z_odd=results_mcmc_odd['50th'][2]
        rp_z_err_odd16=results_mcmc_odd['16th'][2]
        rp_z_err_odd84=results_mcmc_odd['84th'][2]
        
        rp_r_odd_med=results_mcmc_odd['Median'][1]
        rp_z_odd_med=results_mcmc_odd['Median'][2]
    
    rp_r_all_batman=results_batman_all['rp_r']
    rp_r_all_err_batman=results_batman_all['rp_r_err']
    rp_z_all_batman=results_batman_all['rp_z']
    rp_z_all_err_batman=results_batman_all['rp_z_err']
    rp_r_all=results_mcmc_all['50th'][1]
    rp_r_err_all16=results_mcmc_all['16th'][1]
    rp_r_err_all84=results_mcmc_all['84th'][1] 
    rp_z_all=results_mcmc_all['50th'][2]
    rp_z_err_all16=results_mcmc_all['16th'][2]
    rp_z_err_all84=results_mcmc_all['84th'][2]
    
    rp_r_all_med=results_mcmc_all['Median'][1]
    rp_z_all_med=results_mcmc_all['Median'][2]
    
    if results_mcmc_even is not None and results_mcmc_odd is not None:
        x_r = [1, 2, 3]
        x_z = [1.1, 2.1, 3.1]
        y_r = [rp_r_all, rp_r_even, rp_r_odd]
        y_z = [rp_z_all, rp_z_even, rp_z_odd]
        y_r_batman = [rp_r_all_batman, rp_r_even_batman, rp_r_odd_batman]
        y_z_batman = [rp_z_all_batman, rp_z_even_batman, rp_z_odd_batman]
        err_r_batman = [rp_r_all_err_batman, rp_r_even_err_batman, rp_r_odd_err_batman]
        err_z_batman = [rp_z_all_err_batman, rp_z_even_err_batman, rp_z_odd_err_batman]

        err_r_lower = [rp_r_err_all16, rp_r_err_even16, rp_r_err_odd16]
        err_r_upper = [rp_r_err_all84, rp_r_err_even84, rp_r_err_odd84]
        asymmetric_error_r = np.array(list(zip(err_r_lower, err_r_upper))).T

        err_z_lower = [rp_z_err_all16, rp_z_err_even16, rp_z_err_odd16]
        err_z_upper = [rp_z_err_all84, rp_z_err_even84, rp_z_err_odd84]
        asymmetric_error_z = np.array(list(zip(err_z_lower, err_z_upper))).T

        y_r_median = [rp_r_all_med, rp_r_even_med, rp_r_odd_med]
        y_z_median = [rp_z_all_med, rp_z_even_med, rp_z_odd_med]
    elif results_mcmc_even is None and results_mcmc_odd is not None:
        x_r = [1, 3]
        x_z = [1.1, 3.1]
        y_r = [rp_r_all, rp_r_odd]
        y_z = [rp_z_all, rp_z_odd]
        y_r_batman = [rp_r_all_batman, rp_r_odd_batman]
        y_z_batman = [rp_z_all_batman, rp_z_odd_batman]
        err_r_batman = [rp_r_all_err_batman, rp_r_odd_err_batman]
        err_z_batman = [rp_z_all_err_batman, rp_z_odd_err_batman]

        err_r_lower = [rp_r_err_all16, rp_r_err_odd16]
        err_r_upper = [rp_r_err_all84, rp_r_err_odd84]
        asymmetric_error_r = np.array(list(zip(err_r_lower, err_r_upper))).T

        err_z_lower = [rp_z_err_all16, rp_z_err_odd16]
        err_z_upper = [rp_z_err_all84, rp_z_err_odd84]
        asymmetric_error_z = np.array(list(zip(err_z_lower, err_z_upper))).T

        y_r_median = [rp_r_all_med, rp_r_odd_med]
        y_z_median = [rp_z_all_med, rp_z_odd_med]
    elif results_mcmc_odd is None and results_mcmc_even is not None:
        x_r = [1, 2]
        x_z = [1.1, 2.1]
        y_r = [rp_r_all, rp_r_even]
        y_z = [rp_z_all, rp_z_even]
        y_r_batman = [rp_r_all_batman, rp_r_even_batman]
        y_z_batman = [rp_z_all_batman, rp_z_even_batman]
        err_r_batman = [rp_r_all_err_batman, rp_r_even_err_batman]
        err_z_batman = [rp_z_all_err_batman, rp_z_even_err_batman]

        err_r_lower = [rp_r_err_all16, rp_r_err_even16]
        err_r_upper = [rp_r_err_all84, rp_r_err_even84]
        asymmetric_error_r = np.array(list(zip(err_r_lower, err_r_upper))).T

        err_z_lower = [rp_z_err_all16, rp_z_err_even16]
        err_z_upper = [rp_z_err_all84, rp_z_err_even84]
        asymmetric_error_z = np.array(list(zip(err_z_lower, err_z_upper))).T

        y_r_median = [rp_r_all_med, rp_r_even_med]
        y_z_median = [rp_z_all_med, rp_z_even_med]
    
    
    fig = plt.figure()
    ax = fig.add_axes([0.1, 0.1, 0.8, 0.8]) # main axes
    ax.plot(x_r, y_r_median, '*', color='black', label='median values', zorder=3)
    ax.plot(x_z, y_z_median, '*', color='black', zorder=3)
    ax.errorbar(x_r, y_r, asymmetric_error_r, fmt='o', color='blue', label='r-band', zorder=1)
    ax.errorbar(x_z, y_z, asymmetric_error_z, fmt='o', color='red', label='z-band', zorder=1)
    ax.scatter(x_r, y_r_batman,facecolors='none', edgecolors='dimgrey', label='curve fit values', zorder=2)
    ax.scatter(x_z, y_z_batman,facecolors='none', edgecolors='dimgrey', zorder=2) 
    
    if results_mcmc_even is not None and results_mcmc_odd is not None:
        x_even=[4, 5]
        x_odd=[4.1, 5.1]
        y_even=[rp_r_even, rp_z_even]
        y_odd=[rp_r_odd, rp_z_odd]
        y_even_batman=[rp_r_even_batman, rp_z_even_batman]
        y_odd_batman=[rp_r_odd_batman, rp_z_odd_batman]
        err_even_batman=[rp_r_even_err_batman, rp_z_even_err_batman]
        err_odd_batman=[rp_r_odd_err_batman, rp_z_odd_err_batman]

        err_even16=[rp_r_err_even16,rp_z_err_even16]
        err_even84=[rp_r_err_even84,rp_z_err_even84]
        err_odd16=[rp_r_err_odd16, rp_z_err_odd16]
        err_odd84=[rp_r_err_odd84, rp_z_err_odd84]

        y_even_median=[rp_r_even_med, rp_z_even_med]
        y_odd_median=[rp_r_odd_med, rp_z_odd_med]

        ax.plot([4, 5], y_even_median,'*', color='black', zorder=3)
        ax.plot([4.1, 5.1], y_odd_median,'*', color='black', zorder=3)
        ax.errorbar(x_even, y_even, np.array(list(zip(err_even16, err_even84))).T, fmt='o', color='purple', label='even', zorder=1)
        ax.errorbar(x_odd, y_odd, np.array(list(zip(err_odd16, err_odd84))).T, fmt='o', color='green', label='odd', zorder=1)
        ax.scatter([4, 5], y_even_batman, facecolors='none', edgecolors='dimgrey', zorder=2)
        ax.scatter([4.1, 5.1], y_odd_batman, facecolors='none', edgecolors='dimgrey', zorder=2)
    
    ax.set_ylabel('Planet radii (rp/r$_\star$)')
    if period==None:
        ax.set_title('Comparison of planet radii for folded lightcurves')
    else: 
        ax.set_title('Comparison of planet radii for P='+str(round(period,5)))
        #ax.set_title('Comparison of Planet Radii for Candidate 2')

    #ax.set_ylim(0)
    ax.set_xticks([1, 2, 3, 4, 5])
    ax.set_xticklabels(['all transits: r- and z- \ncomparison', '\n\neven transits: r- and z- \ncomparison','odd transits: r- and z- \ncomparison', '\n\nr-band: odd and even \ncomparision', 'z-band: odd and even \ncomparision'])
    #ax.set_xticklabels(['even','odd','total', 'r-band', 'z-band'])
    ax.legend(bbox_to_anchor=(1.05, 1.0), loc='upper left')
    #ax.set_yticks([-1,0,1])
    #fig.savefig(os.path.join('poster_images', 'comparison_of_rp'+str(period)+'.png'),bbox_inches='tight', dpi=300)
    fig.savefig(os.path.join('figs', 'comparison_of_rp.pdf'),bbox_inches='tight')
    

    
def initial_guesses(period, transit_data, data_split):
    t0_guesses=transit_data['tcbest[n]'].tolist()
    rp_guess = np.mean(np.sqrt((-transit_data['best_depth[n]'])))
    #a_guess = ((2/np.sqrt(4*np.pi))*(period/np.mean(transit_data['best_duration[n]'])))**(2/3)
    a_guess = period/(np.pi*np.mean(transit_data['best_duration[n]']))
    #C_guess = np.mean(transit_data['median-all-data'])
    t0_guess = t0_guesses
    rp_r_guess = rp_guess
    rp_z_guess = rp_guess
    a_guess = a_guess
    inc_guess = 89.
    
    all_mag_r_data = [mag for sublist in data_split['mag_r'].tolist() for mag in sublist]
    all_mag_z_data = [mag for sublist in data_split['mag_z'].tolist() for mag in sublist]
    C_r_guess = np.mean(all_mag_r_data)
    C_z_guess = np.mean(all_mag_z_data)
    p0_guess=[t0_guess, rp_r_guess, rp_z_guess, a_guess, inc_guess, C_r_guess, C_z_guess]
    return p0_guess
    
def load_data(filename_data, starid, exclude=None, exclude_transit=None, include_transit=None):
    columns=['starid', 'lcname', 'starflag', 'nightflag', 'int(night)', 'median-all-data', 'std-all-data', 'tcbest[n]',
    'best_depth[n]', 'best_depth_sigma[n]', 'best_duration[n]', 'snrmax[n]', 'bestmedian_intransit[n]',
    'best_std_intransit[n]', 'best_npoints_intransit[n]', 'best_median_outoftransit[n]', 'best_std_outoftransit[n]', 'best_npoints_outoftransit[n]', 'detections_so_far']
    data=pd.read_csv(filename_data, sep='\s+',
                      header=None)
    #print('data')
    #print(starid)
    data.columns=columns
    star_data = data[data['starid']==starid]
    transit_data = star_data[star_data['snrmax[n]']>5.0]
    transit_nights = transit_data['int(night)'].tolist()
    all_nights = star_data['int(night)'].tolist()
    #print('transit nights',transit_nights)
    #print('all nights', all_nights)
    if exclude_transit != None:
        for num in exclude_transit:
            transit_data=transit_data[transit_data['int(night)']!=num]
        transit_nights=[night for night in transit_nights if night not in exclude_transit]
    if exclude != None:
        all_nights = [night for night in all_nights if night not in exclude]
        transit_nights=[night for night in transit_nights if night not in exclude]
    if include_transit != None:
        transit_nights.extend(include_transit)
        transit_nights=sorted(transit_nights)
    transit_data_new=pd.DataFrame({})
    for night in transit_nights:
        df=star_data[star_data['int(night)']==night]
        transit_data_new=pd.concat([transit_data_new, df], ignore_index=True)
    transit_nights=list(set(transit_nights))
    #print('transit data after')
    #print(transit_data_new)
    #print('transit nights after',transit_nights)
    #print('all nights after', all_nights)
    return star_data, transit_data_new, transit_nights, all_nights
    
def weighted_avg_and_std(values, errors):
    """
    Return the weighted average and standard deviation.

    values, weights -- NumPy ndarrays with the same shape.
    """
    weights = 1/np.square(errors)
    average = np.average(values, weights=weights)
    # Fast and numerically precise:
    variance = np.average((values-average)**2, weights=weights)
    return (average, np.sqrt(variance))

def difference(var1, var1_err, var2, var2_err):
    diff = var1-var2
    diff_err = np.sqrt((var1_err**2)+(var2_err**2))
    return diff, diff_err

def odd_even_data(data, t0_all, period):
    t0_all=np.array(t0_all)
    t0_new = t0_all - t0_all[0]
    #print(t0_new)
    div=[round(t/period) for t in t0_new]
    #print(t0_all)
    #print(div)
    
    time_transits_r = data['time_r']
    time_transits_z = data['time_z']
    mag_transits_r = data['mag_r']
    mag_transits_z = data['mag_z']
    err_transits_r = data['err_r']
    err_transits_z = data['err_z']
    
    time_transits_even_r = []
    mag_transits_even_r = []
    err_transits_even_r = []
    time_transits_odd_r = []
    mag_transits_odd_r = []
    err_transits_odd_r = []

    time_transits_even_z = []
    mag_transits_even_z = []
    err_transits_even_z = []
    time_transits_odd_z = []
    mag_transits_odd_z = []
    err_transits_odd_z = []
    
    indices=np.arange(0, len(div))
    for idx in range(len(time_transits_r)):
        if div[idx]%2==0:
            time_transits_even_r.extend(time_transits_r[idx])
            mag_transits_even_r.extend(mag_transits_r[idx])
            err_transits_even_r.extend(err_transits_r[idx])
            time_transits_even_z.extend(time_transits_z[idx])
            mag_transits_even_z.extend(mag_transits_z[idx])
            err_transits_even_z.extend(err_transits_z[idx])
        else:
            time_transits_odd_r.extend(time_transits_r[idx])
            mag_transits_odd_r.extend(mag_transits_r[idx])
            err_transits_odd_r.extend(err_transits_r[idx])
            time_transits_odd_z.extend(time_transits_z[idx])
            mag_transits_odd_z.extend(mag_transits_z[idx])
            err_transits_odd_z.extend(err_transits_z[idx])
            
    data_even_z = pd.DataFrame({'time_z':time_transits_even_z, 'mag_z': mag_transits_even_z, 'err_z': err_transits_even_z})
    data_even_r = pd.DataFrame({'time_r':time_transits_even_r, 'mag_r': mag_transits_even_r, 'err_r': err_transits_even_r})
    data_odd_r = pd.DataFrame({'time_r':time_transits_odd_r, 'mag_r': mag_transits_odd_r, 'err_r': err_transits_odd_r})
    data_odd_z = pd.DataFrame({'time_z':time_transits_odd_z, 'mag_z': mag_transits_odd_z, 'err_z': err_transits_odd_z})
    return data_even_r, data_even_z, data_odd_r, data_odd_z

def get_split_data(filename_r, filename_z): #got rid of transit_night argument
        
    data_z = pd.read_table(filename_z, header=None, delimiter="\s+")
    time_z = data_z[0] 
    mag_z = data_z[1] 
    err_z = data_z[2] 
    time_z=np.array(time_z)
    mag_z=np.array(mag_z)
    err_z=np.array(err_z)
    time_z = time_z - 2458664

    data_r = pd.read_table(filename_r, header=None, delimiter="\s+")
    time_r = data_r[0] 
    mag_r = data_r[1] 
    err_r = data_r[2]
    time_r=np.array(time_r)
    mag_r=np.array(mag_r)
    err_r=np.array(err_r)
    time_r = time_r - 2458664

    time_split_r = []
    mag_split_r = []
    err_split_r = []
    time_split_z = []
    mag_split_z = []
    err_split_z = []
    j=0
    for i in range(0, len(time_r)):
        if i == len(time_r)-1:
            continue
        if (time_r[i+1]-time_r[i])>0.5:
            time_split_r.append(time_r[j+1:i+1])
            mag_split_r.append(mag_r[j+1:i+1])
            err_split_r.append(err_r[j+1:i+1])

            j=i
        if time_r[i+1]==time_r[-1]:
            time_split_r.append(time_r[j+1:])
            mag_split_r.append(mag_r[j+1:])
            err_split_r.append(err_r[j+1:])

    j=0
    for i in range(0, len(time_z)):
        if i == len(time_z)-1:
            continue
        if (time_z[i+1]-time_z[i])>0.5:
            time_split_z.append(time_z[j+1:i+1])
            mag_split_z.append(mag_z[j+1:i+1])
            err_split_z.append(err_z[j+1:i+1])
            j=i
        if time_z[i+1]==time_z[-1]:
            time_split_z.append(time_z[j+1:])
            mag_split_z.append(mag_z[j+1:])
            err_split_z.append(err_z[j+1:])

#     time_by_night_r = {}
#     mag_by_night_r = {}
#     err_by_night_r = {}
#     idx=0
#     for time in time_split_r:
#         night=int(time[0])
#         time_by_night_r[night]=time
#         mag_by_night_r[night]=mag_split_r[idx]
#         err_by_night_r[night]=err_split_r[idx]
#         idx+=1

#     time_by_night_z = {}
#     mag_by_night_z = {}
#     err_by_night_z = {}
#     idx=0
#     for time in time_split_z:
#         night=int(time[0])
#         time_by_night_z[night]=time
#         mag_by_night_z[night]=mag_split_z[idx]
#         err_by_night_z[night]=err_split_z[idx]
#         idx+=1

#     time_transits_r = [time_by_night_r[night] for night in time_by_night_r if night in transit_nights]
#     time_transits_z = [time_by_night_z[night] for night in time_by_night_z if night in transit_nights]
#     mag_transits_r = [mag_by_night_r[night] for night in mag_by_night_r if night in transit_nights]
#     mag_transits_z = [mag_by_night_z[night] for night in mag_by_night_z if night in transit_nights]
#     err_transits_r = [err_by_night_r[night] for night in err_by_night_r if night in transit_nights]
#     err_transits_z = [err_by_night_z[night] for night in err_by_night_z if night in transit_nights]

#     time_transits_r = np.array(time_transits_r, dtype=object)
#     time_transits_z = np.array(time_transits_z, dtype=object)
#     mag_transits_r = np.array(mag_transits_r, dtype=object)
#     mag_transits_z = np.array(mag_transits_z, dtype=object)
#     err_transits_r = np.array(err_transits_r, dtype=object)
#     err_transits_z = np.array(err_transits_z, dtype=object)
    
    #print(len(time_transits_r), len(time_transits_z), len(mag_transits_r), len(mag_transits_z))

    data = pd.DataFrame({'time_r':time_split_r, 'time_z':time_split_z, 'mag_r':mag_split_r, 'mag_z':mag_split_z, 'err_r':err_split_r, 'err_z':err_split_z})
    
    return data

def get_transit_data(data_split, transit_nights):
    #for night in transit_nights:
        #dr = lcr[(lcr.iloc[:,tcol]>night) & (lcr.iloc[:,tcol]<night+1.0)]
        #dz = lcz[(lcz.iloc[:,tcol]>night) & (lcz.iloc[:,tcol]<night+1.0)]
        
    time_split_r=list(data_split['time_r'])
    time_split_z=list(data_split['time_z'])
    mag_split_r=list(data_split['mag_r'])
    mag_split_z=list(data_split['mag_z'])
    err_split_r=list(data_split['err_r'])
    err_split_z=list(data_split['err_z'])
    
    time_by_night_r = {}
    mag_by_night_r = {}
    err_by_night_r = {}
    idx=0
    for time in time_split_r:
        night=int(time[0])
        time_by_night_r[night]=time
        mag_by_night_r[night]=mag_split_r[idx]
        err_by_night_r[night]=err_split_r[idx]
        idx+=1

    time_by_night_z = {}
    mag_by_night_z = {}
    err_by_night_z = {}
    idx=0
    for time in time_split_z:
        night=int(time[0])
        time_by_night_z[night]=time
        mag_by_night_z[night]=mag_split_z[idx]
        err_by_night_z[night]=err_split_z[idx]
        idx+=1

    time_transits_r = [time_by_night_r[night] for night in time_by_night_r if night in transit_nights]
    time_transits_z = [time_by_night_z[night] for night in time_by_night_z if night in transit_nights]
    mag_transits_r = [mag_by_night_r[night] for night in mag_by_night_r if night in transit_nights]
    mag_transits_z = [mag_by_night_z[night] for night in mag_by_night_z if night in transit_nights]
    err_transits_r = [err_by_night_r[night] for night in err_by_night_r if night in transit_nights]
    err_transits_z = [err_by_night_z[night] for night in err_by_night_z if night in transit_nights]

    time_transits_r = np.array(time_transits_r, dtype=object)
    time_transits_z = np.array(time_transits_z, dtype=object)
    mag_transits_r = np.array(mag_transits_r, dtype=object)
    mag_transits_z = np.array(mag_transits_z, dtype=object)
    err_transits_r = np.array(err_transits_r, dtype=object)
    err_transits_z = np.array(err_transits_z, dtype=object)
    
    data = pd.DataFrame({'time_r':time_transits_r, 'time_z':time_transits_z, 'mag_r':mag_transits_r, 'mag_z':mag_transits_z, 'err_r':err_transits_r, 'err_z':err_transits_z})
    
    return data
    
    
def folded_lightcurve(data_r, data_z, period, t0_1, show_plot=True, ans_r=None, ans_z=None, initial_guess_r=None, initial_guess_z=None):
    
    T0 = 0.0
    
    time_r = data_r['time_r']
    mag_r = data_r['mag_r']
    err_r = data_r['err_r']
    time_z = data_z['time_z']
    mag_z = data_z['mag_z']
    err_z = data_z['err_z']
    
    time_r = np.array(time_r) 
    time_z = np.array(time_z)
    
    time_r = time_r - t0_1 + 0.25*period
    time_z = time_z - t0_1 + 0.25*period
    
    phases_r = foldAt(time_r, period, T0)
    sortIndi = np.argsort(phases_r)
    # ... and, second, rearrange the arrays.
    phases_r = phases_r[sortIndi]
    mag_r = np.array(mag_r)[sortIndi]
    err_r = np.array(err_r)[sortIndi]
    if type(ans_r)==np.ndarray:
        ans_r=ans_r[sortIndi]
    if type(initial_guess_r)==np.ndarray:
        initial_guess_r=initial_guess_r[sortIndi]
    
    phases_z = foldAt(np.array(time_z), period, T0)
    sortIndi = np.argsort(phases_z)
    # ... and, second, rearrange the arrays.
    phases_z = phases_z[sortIndi]
    mag_z = np.array(mag_z)[sortIndi]
    err_z = np.array(err_z)[sortIndi]
    if type(ans_z)==np.ndarray:
        ans_z=ans_z[sortIndi]
    if type(initial_guess_z)==np.ndarray:
        initial_guess_z=initial_guess_z[sortIndi]

    if show_plot==True:
        plt.plot(phases_r, mag_r, 'o', color='red', markersize=1.)
        plt.plot(phases_z, mag_z, 'o', color='blue', markersize=1.)
        plt.gca().invert_yaxis()
        plt.title('Folded lightcurve for P='+str(period))
        plt.xlabel('Phase')
        plt.ylabel('Brightness')
        #plt.show()

    data_r_new = pd.DataFrame({'phases_r': phases_r, 'mag_r':mag_r, 'err_r': err_r})
    data_z_new = pd.DataFrame({'phases_z': phases_z, 'mag_z':mag_z, 'err_z': err_z})
    
    if type(ans_r)==np.ndarray and type(ans_z)==np.ndarray:
        return data_r_new, data_z_new, ans_r, ans_z, initial_guess_r, initial_guess_z
    else:
        return data_r_new, data_z_new


# %%
