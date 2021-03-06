/*

 $Id: bestpdb.h,v 1.5 2009/05/08 23:02:10 rhuey Exp $

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

#ifndef BESTPDB
#define BESTPDB

#include "constants.h"
#include "print_rem.h"
#include "strindex.h"
#include "print_avsfld.h"

void  bestpdb( int   ncluster, 
               int   num_in_clu[MAX_RUNS], 
               int   cluster[MAX_RUNS][MAX_RUNS], 
               Real econf[MAX_RUNS], 
               Real crd[MAX_RUNS][MAX_ATOMS][SPACE], 
               char  atomstuff[MAX_ATOMS][MAX_CHARS], 
               int   natom, 
               Boole B_write_all_clusmem, 
               Real ref_rms[MAX_RUNS]);
#endif
