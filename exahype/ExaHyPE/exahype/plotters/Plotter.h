/**
 * This file is part of the ExaHyPE project.
 * Copyright (c) 2016  http://exahype.eu
 * All rights reserved.
 *
 * The project has received funding from the European Union's Horizon 
 * 2020 research and innovation programme under grant agreement
 * No 671698. For copyrights and licensing, please consult the webpage.
 *
 * Released under the BSD 3 Open Source License.
 * For the full license text, see LICENSE.txt
 **/
 
#ifndef _EXAHYPE_PLOTTERS_PLOTTER_H_
#define _EXAHYPE_PLOTTERS_PLOTTER_H_

#include <string>
#include <vector>

#include "exahype/Parser.h"
#include "peano/utils/Globals.h"

#include "tarch/la/Vector.h"
#include "tarch/logging/Log.h"

namespace exahype {
  namespace plotters {
    class Plotter;

    extern std::vector<Plotter*> RegisteredPlotters;

    bool isAPlotterActive(double currentTimeStamp);
    void finishedPlotting();
    double getTimeOfNextPlot();
  }
}


/**
 * Central plotter class
 *
 * The ExaHyPE kernel holds one plotter instance per plotter specified in the
 * config file.
 *
 * @author Tobias Weinzierl
 */
class exahype::plotters::Plotter {
 public:

  /**
   * Interface/abstract superclass for user-defined filtering and in-situ postprocessing
   *
   * For each solver that one specifies in the spec file, the toolkit generates
   * one filter class where the user can in-situ postprocess the data before it
   * goes into the actual output file.
   *
   * The superclass of such an in-situ postprocessing is this type.
   *
   * @author Tobias Weinzierl
   */
  class UserOnTheFlyPostProcessing {
    public:
      virtual ~UserOnTheFlyPostProcessing() {}

      /**
       * Start plotting
       *
       * Is called per plot. The code instantiates each postprocessing filter
       * only once. However, it calls the start and finish routine per written
       * file. This way, you can either keep track of global data over the
       * whole simulation time (if you plug into constructors and desctructors)
       * or you can keep track of global quantities per written snapshot.
       *
       * Time is the minimal global time stamp of the simulation when the plot
       * is started. If you have global time stepping, time is the simulation
       * time of all data. If you use local or anarchic time stepping, time is
       * the minimum of all time steps over the domain, as a plotter becomes
       * active every time the minimum of all data in the domain overruns the
       * snapshot time stamp.
       */
      virtual void startPlotting( double time ) = 0;

      /**
       * Counterpart of startPlotting
       */
      virtual void finishPlotting() = 0;

      /**
       * Mapping of simulation quantities onto output quantities
       *
       * This routine is called per output quantity that is to be plotted.
       * Please note that it is called only for quantities of cells that do
       * overlap with the spatial filter rules (if there are any). Filtering
       * is done on a per-cell (ADER-DG)/per-patch (Finite Volumes) basis, i.e.
       * a quantity might be outside of the filter and be plotted nevertheless
       * if the corresponding cell overlaps with the plotted region.
       *
       * The sample locations (where the quantity comes from) depends on the
       * type of the plotter underlying the filter. If it is a Cartesian
       * plotter, then we sample all higher order polynomials within each cell
       * with regular spacing equal to the order. If you use a probe, the
       * quantity stems exactly from the probe (cmp. x parameter). Consult
       * the spec of the plotters for details if in doubt.
       *
       * @param offsetOfPatch Offset (left, bottom, front corner) of patch from
       *          which the quantity stems from.
       * @param sizeOfPatch   Size of the underlying patch.
       * @param x             Location of the mapped quantity. x is always
       *          contained within the cell specified via offsetOfPatch and
       *          sizeOfPatch. If you use a Cartesian plotter, the x are all
       *          equidistant. Anyway, it is always the exact location where
       *          Q is taken from the solution.
       * @param pos           All plotters break down each ADER-DG cell into
       *          grids. For Finite Volumes, such grids are naturally in place
       *          already. pos it the internal counter without this grids. For
       *          ADER-DG, it enumerates through the Gauss-Legendre nodes. So
       *          in this case, pos and x are slightly redundant - you can
       *          compute x from pos. However, we decided to provide both.
       *          If you invoke a probe plotter, this entry always is the zero
       *          vector.
       * @param Q             Vector of unknowns at position x. The cardinality
       *          of the data is the unknowns plus the material parameters from
       *          the spec file's solver.
       * @param outputQuantities Vector of unknowns to be plotted actually. The
       *          cardinality of this vector is the size specified via unknowns
       *          in the plotter section in the spec file. Mappings are
       *          expected to befill this array with data (usually from Q). If
       *          the output cardinality in the spec file equals 0, this pointer
       *          is nullptr.
       * @param timeStamp     Actual time stamp of the quantity. In local time
       *          stepping, different patches might have different time stamps
       *          when we plot.
       */
      virtual void mapQuantities(
        const tarch::la::Vector<DIMENSIONS, double>& offsetOfPatch,
        const tarch::la::Vector<DIMENSIONS, double>& sizeOfPatch,
        const tarch::la::Vector<DIMENSIONS, double>& x,
	const tarch::la::Vector<DIMENSIONS, int>&    pos,
        double* Q,
        double* outputQuantities,
        double timeStamp) = 0;
  };

  /**
   * Actual device chosen by a plotter in the config file. If you implement
   * your own device, please also add a
   *
   * static std::string getIdentifier();
   *
   * routine that returns the device's identifier.
   */
  class Device {
   protected:
    UserOnTheFlyPostProcessing*  _postProcessing;
   public:
    Device(UserOnTheFlyPostProcessing* postProcessing):
      _postProcessing(postProcessing) {}
    virtual ~Device() {}

    /**
     * Configure the plotter device. Is invoked directly after the constructor is
     * called.
     */
    virtual void init(const std::string& filename, int order, int unknowns, int writtenUnknowns, const std::string& select) = 0;


    /**
     * Hand a patch over to the plotter device.
     */
    virtual void plotPatch(
        const int cellDescriptionsIndex,
        const int element) = 0;

    virtual void startPlotting( double time ) = 0;
    virtual void finishPlotting() = 0;
  };

 private:
  static tarch::logging::Log _log;

  int              _solver;
  std::string      _identifier;
  int              _writtenUnknowns;
  /**
   * Time of the next snapshot to be written
   */
  double                 _time;
  /**
   * Most recent time stamp communicated by the solver.
   */
  double                 _solverTimeStamp;
  const double           _repeat;
  std::string            _filename;
  const std::string      _select;
  bool                   _isActive;

  Device*                _device;

 public:
  /**
   * @param solverConfig Number of the underlying solver. This number is important to
   *               parse the file: the constructor asks parser for the \p solverNumber-th
   *               solver section.
   * @param plotterConfig Same story: Required to tell the parser which tag in
   *                      the file is to be read in.
   */
  Plotter(const int solverConfig,const int plotterConfig,
          const exahype::Parser& parser, UserOnTheFlyPostProcessing* postProcessing);

  /**
   * Plotter constructor for scenarios where we want to use a plotter configuration
   * given in the specification file for another solver than the one
   * it is specified for.
   *
   * @param solverConfig Number of the underlying solver. This number is important to
   *               parse the file: the constructor asks parser for the \p solverNumber-th
   *               solver section.
   * @param plotterConfig Same story: Required to tell the parser which tag in
   *                      the file is to be read in.
   *
   * @param solverDataSource Number of the solver this plotter gets its data from.
   *               For parameter studies we usually want to use the same
   *               plotter configuration for multiple solvers.
   */
  Plotter(const int solverConfig,const int plotterConfig,
          const exahype::Parser& parser, UserOnTheFlyPostProcessing* postProcessing,
          const int solverDataSource);
  ~Plotter();

  // Disallow copy and assignment
  Plotter(const Plotter& other) = delete;
  Plotter& operator=(const Plotter& other) = delete;

  /**
   * Checks whether there should be a plotter according to this class.
   * If it should become open, it is opened
   */
  bool checkWetherPlotterBecomesActive(double currentTimeStamp);
  bool isActive() const;
  bool plotDataFromSolver(int solver) const;

  /**
   * <h2>Large solver time steps</h2>
   * In case of large solver time step sizes, it can
   * happen that two or more plotting times fall into the
   * interval between the current solver time stamp
   * and the previous one.
   *
   * In these scenarios, we have to increment the
   * plotting time multiple times by the
   * repeat time.
   */
  void finishedPlotting();

  void plotPatch(
      const int cellDescriptionsIndex,const int element);

  std::string getFileName() const;

  /**
   * If you want to check whether to active a plotter, please use
   * checkWetherSolverBecomesActive(). This operation is meant to determine
   * the next time step if you have fixed time stepping.
   */
  double getNextPlotTime() const;

  std::string toString() const;
};

#endif
