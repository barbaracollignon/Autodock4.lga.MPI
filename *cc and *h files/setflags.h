/*

 $Id: setflags.h,v 1.6 2009/05/08 23:02:17 rhuey Exp $
 $Id: setflags.h, autodock4.lga.MPI v0.0 2010/10/10 22:55:01 collignon Exp $

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

#ifndef SETFLAGS
#define SETFLAGS
#include "constants.h"
#include "usage.h"

int setflags( int argc, 
              char ** argv,
              char * version_num, 
              int id,
              int id2);

#endif
