/**
 *  $File: config.h $
 *  $LastChangedDate: 2012-03-23 21:40:11 -0500 (Fri, 23 Mar 2012) $
 *  $Rev: 1013 $
 *
 *  This file is part of simuPOP, a forward-time population genetics
 *  simulation environment. Please visit http://simupop.sourceforge.net
 *  for details.
 *
 *  Copyright (C) 2004 - 2010 Bo Peng (bpeng@mdanderson.org)
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program. If not, see <http://www.gnu.org/licenses/>.
 */

/** \file config.h
 *  This file includes appropriate configuration file for different
 *  operating systems and compilers.
 */

#ifdef _WIN32
#  include "config_win32.h"
#else
#  ifdef MACOSX
#    include "config_mac.h"
#  else
#    ifdef __sparc__
#      include "config_solaris.h"
#    else
#      include "config_linux.h"
#    endif
#  endif
#endif
