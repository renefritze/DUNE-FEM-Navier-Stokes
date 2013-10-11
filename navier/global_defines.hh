#ifndef DUNE_NAVIERSTOKES_GLOBAL_DEFINES_HH
#define DUNE_NAVIERSTOKES_GLOBAL_DEFINES_HH

// gcc header cyclic dependency error workaround
namespace std {
class type_info;
}

#define TESTCASE_NAME "TESTCASE"

#if !defined(TESTCASE)
#define TESTCASE TestCase3D
#endif

//#if ( ( defined(SGRID) || defined(ALUGRID_SIMPLEX) ||  defined(ALUGRID_CUBE) ) && ( GRIDDIM == 3 ) ) ||
//defined(UGGRID) || defined(YASPGRID)
//	//this is no mistake, ALU is indeed only incompatible in 3d
//	#define OLD_DUNE_GRID_VERSION
//#endif

#define MODEL_PROVIDES_LOCALFUNCTION 0

#ifndef NAVIER_DATA_NAMESPACE
#define NAVIER_DATA_NAMESPACE NavierProblems::Runtime
#endif

// the adaption manager might be troublesome with certain gridparts/spaces, so we needed a easy way to disable it
#ifndef ENABLE_ADAPTIVE
#define ENABLE_ADAPTIVE 1
#endif

#ifndef COMMIT
#define COMMIT "undefined"
#endif

#include <string>
static const std::string commit_string(COMMIT);

#endif // DUNE_NAVIERSTOKES_GLOBAL_DEFINES_HH

/** Copyright (c) 2012, Rene Milk
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice, this
 *    list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions and the following disclaimer in the documentation
 *    and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
 * ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * The views and conclusions contained in the software and documentation are those
 * of the authors and should not be interpreted as representing official policies,
 * either expressed or implied, of the FreeBSD Project.
**/
