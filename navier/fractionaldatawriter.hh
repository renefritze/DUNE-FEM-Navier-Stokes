#ifndef FRACTIONALDATAWRITER_HH
#define FRACTIONALDATAWRITER_HH

#include <dune/fem/io/file/datawriter.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <dune/stuff/magnitudefunction.hh>
#include "fractionaltimeprovider.hh"

namespace Dune {
namespace NavierStokes {
/** \brief a Dune::DataWriter derivative that handles current time internaly via passed Dune::Timeprovider
    \note TimeproviderType can only be FractionalTimeProvider
    **/
template <class TimeproviderType, class GridType, class OutputTupleType>
class TimeAwareDataWriter : public DataWriter<GridType, OutputTupleType> {
protected:
  typedef DataWriter<GridType, OutputTupleType> BaseType;
  using BaseType::grid_;
  using BaseType::path_;
  using BaseType::datapref_;
  using BaseType::writeStep_;
  using BaseType::outputFormat_;
  using BaseType::vtk;
  using BaseType::vtkvtx;
  using BaseType::gnuplot;
  using BaseType::data_;
  using BaseType::saveTime_;
  using BaseType::saveStep_;

public:
  TimeAwareDataWriter(const TimeproviderType& timeprovider, const GridType& grid, OutputTupleType& tuple)
    : BaseType(grid, tuple, timeprovider)
    , timeprovider_(timeprovider) {}

  void write() const { write(timeprovider_.time(), timeprovider_.timeStep()); }

  /** \brief write given data to disc
     \param[in] time actual time of computation
     \param[in] timestep current number of time step
     \param[in] data data to write (template type, should be a tuple)
  */
  void write(double time, int /*timestep*/) const {
    //				if( BaseType::willWrite( time, timestep ) )
    {
      {
        // check online display
        BaseType::display();
      }

      if (outputFormat_ == vtk || outputFormat_ == vtkvtx) {
        // write data in vtk output format
        writeVTKOutput(get<0>(data_));
      } else if (outputFormat_ == gnuplot) {
        writeGnuPlotOutput(time);
      } else {
        DUNE_THROW(NotImplemented, "DataWriter::write: wrong output format");
      }

      // only write info for proc 0, otherwise on large number of procs
      // this is to much output
      //					if(myRank_ <= 0)
      //					{
      //						std::cout << "DataWriter["<<myRank_<<"]::write:  time = "<< time << "  writeData: step number " <<
      //writeStep_ << std::endl;
      //					}
      saveTime_ += saveStep_;
      ++writeStep_;
    }
    return;
  }

protected:
  const TimeproviderType& timeprovider_;
  static inline std::string genFilename(const std::string& path, const std::string& fn, int ntime, int precision = 6) {
    std::ostringstream name;

    if (path.size() > 0) {
      name << path;
      name << "/";
    }
    name << fn;

    char cp[256];
    switch (precision) {
      case 2: {
        sprintf(cp, "%02d", ntime);
        break;
      }
      case 3: {
        sprintf(cp, "%03d", ntime);
        break;
      }
      case 4: {
        sprintf(cp, "%04d", ntime);
        break;
      }
      case 5: {
        sprintf(cp, "%05d", ntime);
        break;
      }
      case 6: {
        sprintf(cp, "%06d", ntime);
        break;
      }
      case 7: {
        sprintf(cp, "%07d", ntime);
        break;
      }
      case 8: {
        sprintf(cp, "%08d", ntime);
        break;
      }
      case 9: {
        sprintf(cp, "%09d", ntime);
        break;
      }
      case 10: {
        sprintf(cp, "%010d", ntime);
        break;
      }
      default: { DUNE_THROW(Exception, "Couldn't gernerate filename with precision = " << precision); }
    }
    name << cp;

    // here implicitly a string is generated
    return name.str();
  }

  template <class VTKOut>
  class VTKOutputter {
  public:
    //! Constructor
    VTKOutputter(VTKOut& vtkOut, std::string path, bool parallel, int step)
      : vtkOut_(vtkOut)
      , path_(path)
      , parallel_(parallel)
      , step_(step) {}

    //! Applies the setting on every DiscreteFunction/LocalFunction pair.
    template <class DFType>
    void visit(DFType* f) {
      if (!f)
        return;

      // needs to in same scope as clear() ?
      Stuff::MagnitudeFunction<DFType> magnitude(*f);

      std::string name = genFilename((parallel_) ? "" : path_, f->name(), step_);
      if (DFType::FunctionSpaceType::DimRange > 1) {
        vtkOut_.addVectorVertexData(*f);
        vtkOut_.addVectorCellData(*f);
        vtkOut_.addVertexData(magnitude.discreteFunction());
        vtkOut_.addCellData(magnitude.discreteFunction());
      } else {
        vtkOut_.addVertexData(*f);
        vtkOut_.addCellData(*f);
      }
      const bool binary_output = Parameters().getParam("binary_vtk", true);

      if (parallel_) {
        // write all data for parallel runs
        vtkOut_.pwrite(name.c_str(), path_.c_str(), ".", binary_output ? Dune::VTK::base64 : Dune::VTK::ascii);
      } else {
        // write all data serial
        vtkOut_.write(name.c_str(), binary_output ? Dune::VTK::base64 : Dune::VTK::ascii);
      }

      vtkOut_.clear();
    }

  private:
    VTKOut& vtkOut_;
    const std::string path_;
    const bool parallel_;
    const int step_;
  };

  template <class DFType>
  void writeVTKOutput(const DFType* func) const {
    if (!func)
      return;

    // check whether we have parallel run
    const bool parallel = (grid_.comm().size() > 1);

    // generate filename, with path only for serial run
    {
      // get grid part
      typedef typename DFType::DiscreteFunctionSpaceType::GridPartType GridPartType;
      const GridPartType& gridPart = func->space().gridPart();

      {
        // create vtk output handler
        typedef SubsamplingVTKIO<GridPartType> VTKIOType;
        VTKIOType vtkio(gridPart, VTKOptions::nonconforming);

        // add all functions
        ForEachValue<OutputTupleType> forEach(data_);
        VTKOutputter<VTKIOType> io(vtkio, path_, parallel, timeprovider_.timeStep());
        forEach.apply(io);

        // write all data
      }
    }
  }

  struct Gnu {
    const double time_;
    const std::string path_, datapref_;
    const bool parallel_;
    const int step_;
    Gnu(const double time, std::string path, bool parallel, int step, std::string datapref)
      : time_(time)
      , path_(path)
      , datapref_(datapref)
      , parallel_(parallel)
      , step_(step) {}
    // write to gnuplot file format
    template <class DFType>
    void visit(const DFType* func) const {
      if (!func)
        return;

      typedef typename DFType::Traits Traits;
      typedef typename Traits::LocalFunctionType LocalFunctionType;
      typedef typename Traits::DiscreteFunctionSpaceType DiscreteFunctionSpaceType;
      typedef typename DiscreteFunctionSpaceType::IteratorType IteratorType;
      typedef typename DiscreteFunctionSpaceType::GridPartType GridPartType;

      typedef typename DiscreteFunctionSpaceType::DomainType DomainType;
      typedef typename DiscreteFunctionSpaceType::RangeType RangeType;

      enum {
        dimDomain = DiscreteFunctionSpaceType::dimDomain
      };
      enum {
        dimRange = DiscreteFunctionSpaceType::dimRange
      };

      // generate filename
      //					std::string name = genFilename( path_, datapref_, step_ );
      std::string name = genFilename(path_, datapref_, 0);
      name += "_" + func->name() + ".gnu";
      const bool first = (time_ > 0.0);
      std::ios_base::openmode mode = first ? std::ios_base::app : std::ios_base::out;
      std::ofstream gnuout(name.c_str(), mode);

      if (first) {
        gnuout << "#";
        for (int i = 0; i < dimDomain; ++i)
          gnuout << "x_" << i << " ";
        for (int i = 0; i < dimRange; ++i)
          gnuout << "f_" << i << " ";
        gnuout << "time"
               << "\n";
      }
      // start iteration
      IteratorType endit = func->space().end();
      for (IteratorType it = func->space().begin(); it != endit; ++it) {
        CachingQuadrature<GridPartType, 0> quad(*it, func->space().order());
        LocalFunctionType lf = func->localFunction(*it);
        for (size_t i = 0; i < quad.nop(); ++i) {
          RangeType u;
          DomainType x = it->geometry().global(quad.point(i));
          lf.evaluate(quad[i], u);
          for (int i = 0; i < dimDomain; ++i)
            gnuout << x[i] << " ";
          for (int i = 0; i < dimRange; ++i)
            gnuout << u[i] << " ";
          gnuout << time_ << "\n";
        }
      }
      gnuout << "\n\n";
    }
  };

  void writeGnuPlotOutput(const double time) const {
    const bool parallel = (grid_.comm().size() > 1);
    ForEachValue<OutputTupleType> forEach(data_);
    Gnu io(time, path_, parallel, timeprovider_.timeStep(), datapref_);
    forEach.apply(io);
  }
};
} // end namespace NavierStokes
} // end namespace Dune

#endif // FRACTIONALDATAWRITER_HH

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
