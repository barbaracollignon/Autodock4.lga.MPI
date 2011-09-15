/*

 $Id: parse_dpf_line.cc,v 1.23 2009/05/08 23:02:15 rhuey Exp $

 AutoDock 

Copyright (C) 2009 The Scripps Research Institute. All rights reserved.

 AutoDock is a Trade Mark of The Scripps Research Institute.

 This program is free software; you can redistribute it and/or
 modify it under the terms of the GNU General Public License
 as published by the Free Software Foundation; either version 2
 of the License, or (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program; if not, write to the Free Software
 Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include "parse_dpf_line.h"



int parse_dpf_line( char line[LINE_LEN] )

/******************************************************************************/
/*      Name: parse_dpf_line                                                  */
/*  Function: Parse the docking parameter file line                           */
/*Copyright (C) 2009 The Scripps Research Institute. All rights reserved. */
/*----------------------------------------------------------------------------*/
/*    Author: Garrett Morris, The Scripps Research Institute                  */
/*      Date: 19/05/94                                                        */
/*----------------------------------------------------------------------------*/
/*    Inputs: line                                                            */
/*   Returns: integer token describing the keyword found.                     */
/*   Globals: none.                                                           */
/*----------------------------------------------------------------------------*/
/* Modification Record                                                        */
/* Date     Inits   Comments                                                  */
/* 06/09/95 RSH     Changed to an array implementation                        */
/* 19/05/94 GMM     Entered code.                                             */
/******************************************************************************/

{
    int j, i, token = DPF_;               /* return -1 if nothing is recognized. */
    char c[LINE_LEN];

    const struct {
       char *lexeme;
       int tokenvalue;
    } tokentable[] = {{"read_target",DPF_READ_TARGET},
                      {"read_ligand_list", DPF_READ_LIGAND},
                      {"ligand", DPF_MOVE},  
                      {"fld", DPF_FLD}, 
                      {"map", DPF_MAP}, 
                      {"move", DPF_MOVE}, 
                      {"about", DPF_ABOUT}, 
                      {"tran0", DPF_TRAN0}, 
                      {"quat0", DPF_QUAT0}, 
                      {"ndihe", DPF_NDIHE}, 
                      {"dihe0", DPF_DIHE0}, 
                      {"torsdof", DPF_TORSDOF}, 
                      {"tstep", DPF_TSTEP}, 
                      {"qstep", DPF_QSTEP}, 
                      {"dstep", DPF_DSTEP}, 
                      {"trnrf", DPF_TRNRF}, 
                      {"quarf", DPF_QUARF}, 
                      {"dihrf", DPF_DIHRF}, 
                      {"flex", DPF_FLEX}, 
                      {"intnbp_coeffs", DPF_INTNBP_COEFFS}, 
                      {"rt0", DPF_RT0}, 
                      {"rtrf", DPF_RTRF}, 
                      {"runs", DPF_RUNS}, 
                      {"cycles", DPF_CYCLES}, 
                      {"accs", DPF_ACCS}, 
                      {"rejs", DPF_REJS}, 
                      {"select", DPF_SELECT}, 
                      {"outlev", DPF_OUTLEV}, 
                      {"rmstol", DPF_RMSTOL}, 
                      {"trjfrq", DPF_TRJFRQ}, 
                      {"trjbeg", DPF_TRJBEG}, 
                      {"trjend", DPF_TRJEND}, 
                      {"trjout", DPF_TRJOUT}, 
                      {"trjsel", DPF_TRJSEL}, 
                      {"extnrg", DPF_EXTNRG}, 
                      {"newcrd", DPF_NEWCRD}, 
                      {"cluster", DPF_CLUSTER}, 
                      {"write_all", DPF_CLUSALL}, 
                      {"write_all_cluster_members", DPF_CLUSALL}, 
                      {"charmap", DPF_CHARMAP}, 
                      {"rmsnosym", DPF_RMSNOSYM}, 
                      {"rmsref", DPF_RMSREF}, 
                      {"watch", DPF_WATCH}, 
                      {"linear_schedule", DPF_SCHEDLIN}, 
                      {"schedule_linear", DPF_SCHEDLIN}, 
                      {"linsched", DPF_SCHEDLIN}, 
                      {"schedlin", DPF_SCHEDLIN}, 
                      {"intelec", DPF_INTELEC}, 
                      {"seed", DPF_SEED}, 
                      {"e0max", DPF_E0MAX}, 
                      {"simanneal", DPF_SIMANNEAL}, 
                      {"hardtorcon", DPF_HARDTORCON}, 
                      {"intnbp_r_eps", DPF_INTNBP_REQM_EPS}, 
                      {"gausstorcon", DPF_GAUSSTORCON}, 
                      {"barrier", DPF_BARRIER}, 
                      {"showtorpen", DPF_SHOWTORPEN}, 
                      {"ga_run", DPF_GALS}, 
                      {"gals_run", DPF_GALS}, 
                      {"do_gals", DPF_GALS}, 
                      {"set_ga", DPF_SET_GA}, 
                      {"set_sw1", DPF_SET_SW1}, 
                      {"set_psw1", DPF_SET_PSW1}, 
                      {"analysis", DPF_ANALYSIS}, 
                      {"ga_pop_size", GA_pop_size}, 
                      {"ga_num_generations", GA_num_generations}, 
                      {"ga_num_evals", GA_num_evals}, 
                      {"ga_window_size", GA_window_size}, 
                      {"ga_low", GA_low}, 
                      {"ga_high", GA_high}, 
                      {"ga_elitism", GA_elitism}, 
                      {"ga_mutation_rate", GA_mutation_rate}, 
                      {"ga_crossover_rate", GA_crossover_rate}, 
                      {"ga_cauchy_alpha", GA_Cauchy_alpha}, 
                      {"ga_cauchy_beta", GA_Cauchy_beta}, 
                      {"sw_max_its", SW_max_its}, 
                      {"sw_max_succ", SW_max_succ}, 
                      {"sw_max_fail", SW_max_fail}, 
                      {"sw_rho", SW_rho}, 
                      {"sw_lb_rho", SW_lb_rho}, 
                      {"psw_trans_scale", PSW_TRANS_SCALE}, 
                      {"psw_rot_scale", PSW_ROT_SCALE}, 
                      {"psw_tors_scale", PSW_TORS_SCALE}, 
                      {"do_local_only", DPF_LS}, 
                      {"ls_run", DPF_LS}, 
                      {"do_global_only", DPF_GS}, 
                      {"ga_only_run", DPF_GS}, 
                      {"ls_search_freq", LS_search_freq}, 
                      {"bin_energies_by_rmsd", DPF_INVESTIGATE}, 
                      {"investigate", DPF_INVESTIGATE}, 
              {"ligand_is_not_inhibitor", DPF_LIG_NOT_INHIB}, 
              {"template", DPF_TEMPL_ENERGY}, 
              {"template_energy_file", DPF_TEMPL_ENERGY}, 
              {"include_1_4_interactions", DPF_INCLUDE_1_4_INTERACTIONS}, 
              {"parameter_library", DPF_PARAMETER_LIBRARY}, 
              {"parameter_file", DPF_PARAMETER_LIBRARY} 
              , {"receptor_types", DPF_RECEPTOR_TYPES}  
              , {"ligand_types", DPF_LIGAND_TYPES}      
              , {"unbound", DPF_UNBOUND}      
              , {"epdb", DPF_EPDB}      
              , {"ga_termination_criterion", DPF_TERMINATION}      
              , {"ga_termination", DPF_TERMINATION}      
              , {"ga_crossover_mode", GA_CROSSOVER_MODE}      
              , {"output_pop_file", DPF_POPFILE}      
              , {"set_pattern", DPF_SET_PATTERN}      
              , {"compute_unbound_extended", DPF_COMPUTE_UNBOUND_EXTENDED} 
              , {"set_unbound_energy", DPF_UNBOUND}      
              , {"flexible_residues", DPF_FLEXRES} 
              , {"flexres", DPF_FLEXRES} 
              , {"elecmap", DPF_ELECMAP} 
              , {"desolvmap", DPF_DESOLVMAP} 
              , {"unbound_intnbp_coeffs", DPF_UNBOUND_INTNBP_COEFFS} 
              , {"rmsatoms", DPF_RMSATOMS} 
              , {"confsampler", DPF_CONFSAMPLER} 
              , {"reorient", DPF_REORIENT} 
              , {"axisangle0", DPF_AXISANGLE0} 
              , {"quaternion0", DPF_QUATERNION0} 
              , {"copyright", DPF_COPYRIGHT} 
              , {"warranty", DPF_WARRANTY} 
              , {"autodock_parameter_version", DPF_PARAMETER_VERSION} 
              , {"unbound_model", DPF_UNBOUND_MODEL} 
              , {"unbound_energy", DPF_UNBOUND} 

#if defined(USING_COLINY)
              , {"coliny", DPF_COLINY}  
#endif
              , {"", DPF_NULL}
              };

    // build lower-case version of the token into 'c'
    c[0] = '\0';
    for (j=0; line[j]!='\0' && !isspace(line[j]); j++) {
        /*  Ignore case */
        c[j] = (char)tolower((int)line[j]);
        /*(void)fprintf(stderr,"%c",c[j]);*/
    }
    /*(void)fprintf(stderr,"/n,j = %d\n",j);*/

    /*  Recognize one character tokens  */

    if ((c[0]=='\n') || (c[0]=='\0')) {
        token = DPF_NULL;
    } else if (c[0]=='#') {
        token = DPF_COMMENT;
    } else for (i=0;  (tokentable[i].tokenvalue!=DPF_NULL) ;  i++) {
    /*  Recognize token strings  */
        /*(void)fprintf(stderr,"i = %d, tokentable[i].lexeme = %s, tokentable[i].value = %d, c = %s\n",i,tokentable[i].lexeme,tokentable[i].tokenvalue,c);*/
        if (strncasecmp(tokentable[i].lexeme, c, j) == 0) {
            token = tokentable[i].tokenvalue;
            break;
        }
    }
    return(token);
}
/* EOF */
