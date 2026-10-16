/* Copyright 2017-2021 Institute for Automation of Complex Power Systems,
 *                     EONERC, RWTH Aachen University
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 *********************************************************************************/

#pragma once

#include <dpsim-models/Base/Base_Ph1_PiLine.h>
#include <dpsim-models/CompositePowerComp.h>
#include <dpsim-models/SP/SP_Ph1_Capacitor.h>
#include <dpsim-models/SP/SP_Ph1_Inductor.h>
#include <dpsim-models/SP/SP_Ph1_Resistor.h>
#include <dpsim-models/Solver/MNATearInterface.h>
#include <dpsim-models/Solver/PFSolverInterfaceBranch.h>

namespace CPS {
namespace SP {
namespace Ph1 {

/// \brief PI-line static phasor model
///
/// For MNA, the series branch is represented according to the actual
/// parameters:
///   R > 0, L > 0 : Resistor -- virtual node -- Inductor
///   R = 0, L > 0 : Inductor directly between the terminals
///   R > 0, L = 0 : Resistor directly between the terminals
///   R = 0, L = 0 : unsupported ideal short; rejected explicitly
///
/// This avoids creating zero-valued series subcomponents whose admittance
/// would otherwise become infinite.
class PiLine : public CompositePowerComp<Complex>,
               public Base::Ph1::PiLine,
               public MNATearInterface,
               public SharedFactory<PiLine>,
               public PFSolverInterfaceBranch {
private:
  /// base voltage [V]
  Real mBaseVoltage;

  /// true if mParallelCap was never given a real value in setParameters()
  /// and holds the MNA-conditioning fallback instead
  bool mParallelCapIsFallback = false;

  /// true if mParallelCond was never given a real value in setParameters()
  /// and holds the MNA-conditioning fallback instead
  bool mParallelCondIsFallback = false;

public:
  // #### Power flow results ####

  /// branch current flow [A], coeff(0) from node 0, coeff(1) from node 1
  const Attribute<MatrixComp>::Ptr mCurrent;

  /// branch active power flow [W]
  const Attribute<Matrix>::Ptr mActivePowerBranch;

  /// branch reactive power flow [Var]
  const Attribute<Matrix>::Ptr mReactivePowerBranch;

protected:
  /// base current [A]
  Real mBaseCurrent;

  /// base apparent power [VA]
  Real mBaseApparentPower;

  /// base omega [1/s]
  Real mBaseOmega;

  /// base impedance [Ohm]
  Real mBaseImpedance;

  /// base admittance [S]
  Real mBaseAdmittance;

  /// base inductance [H]
  Real mBaseInductance;

  /// base capacitance [F]
  Real mBaseCapacitance;

  /// resistance [pu]
  Real mSeriesResPerUnit;

  /// capacitance [pu]
  Real mParallelCapPerUnit;

  /// inductance [pu]
  Real mSeriesIndPerUnit;

  /// conductance [pu]
  Real mParallelCondPerUnit;

  // #### Admittance matrix stamp ####

  MatrixComp mY_element;

  /// Optional series inductance submodel
  std::shared_ptr<Inductor> mSubSeriesInductor;

  /// Optional series resistance submodel
  std::shared_ptr<Resistor> mSubSeriesResistor;

  /// Parallel resistor submodel at terminal 0
  std::shared_ptr<Resistor> mSubParallelResistor0;

  /// Parallel capacitor submodel at terminal 0
  std::shared_ptr<Capacitor> mSubParallelCapacitor0;

  /// Parallel resistor submodel at terminal 1
  std::shared_ptr<Resistor> mSubParallelResistor1;

  /// Parallel capacitor submodel at terminal 1
  std::shared_ptr<Capacitor> mSubParallelCapacitor1;

  /// Right side vectors of subcomponents
  std::vector<const Matrix *> mRightVectorStamps;

public:
  /// nodal active power injection
  const Attribute<Real>::Ptr mActivePowerInjection;

  /// nodal reactive power injection
  const Attribute<Real>::Ptr mReactivePowerInjection;

  // #### General ####

  /// Defines UID, name and logging level
  PiLine(String uid, String name, Logger::Level logLevel = Logger::Level::off);

  /// Defines name and logging level
  PiLine(String name, Logger::Level logLevel = Logger::Level::off)
      : PiLine(name, name, logLevel) {}

  SimPowerComp<Complex>::Ptr clone(String copySuffix) override;

  void setParameters(Real resistance, Real inductance, Real capacitance = -1,
                     Real conductance = -1);

  /// Constructs and registers MNA subcomponents; idempotent.
  void createSubComponents() override;

  /// Derives values from power flow data and initializes internal values.
  void initializeParentFromNodesAndTerminals(Real frequency) override;

  // #### Powerflow section ####

  Real getBaseVoltage() const;

  void setBaseVoltage(Real baseVoltage);

  void calculatePerUnitParameters(Real baseApparentPower, Real baseOmega);

  void pfApplyAdmittanceMatrixStamp(SparseMatrixCompRow &Y) override;

  void updateBranchFlow(VectorComp &current, VectorComp &powerflow);

  void storeNodalInjection(Complex powerInjection);

  // #### Getter ####

  MatrixComp Y_element();

  // #### MNA section ####

  void mnaCompUpdateCurrent(const Matrix &leftVector) override;

  void mnaCompUpdateVoltage(const Matrix &leftVector) override;

  void mnaParentPostStep(Real time, Int timeStepCount,
                         Attribute<Matrix>::Ptr &leftVector) override;

  void
  mnaParentAddPostStepDependencies(AttributeBase::List &prevStepDependencies,
                                   AttributeBase::List &attributeDependencies,
                                   AttributeBase::List &modifiedAttributes,
                                   Attribute<Matrix>::Ptr &leftVector) override;

  // #### MNA tear section ####

  MNAInterface::List mnaTearGroundComponents() override;

  void mnaTearInitialize(Real omega, Real timeStep) override;

  void mnaTearApplyMatrixStamp(SparseMatrixRow &tearMatrix) override;

  void mnaTearApplyVoltageStamp(Matrix &voltageVector) override;

  void mnaTearPostStep(Complex voltage, Complex current) override;
};

} // namespace Ph1
} // namespace SP
} // namespace CPS
