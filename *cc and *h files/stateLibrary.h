/*

 $Id: stateLibrary.h,v 1.5 2009/05/08 23:02:17 rhuey Exp $

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

#ifndef COPYSTATE
#define COPYSTATE

#include "structs.h"
#include "constants.h"

void initialiseState( State *S );

void initialiseQuat( Quat *Q );

void copyState( State *destination,
		State  source);

void printState( FILE *fp,
		 State state, 
		 int detail );

void writeState( FILE *fp, 
		 State state );

int checkState( const State *D );

Molecule copyStateToMolecule(State *source, Molecule *mol);
#endif
