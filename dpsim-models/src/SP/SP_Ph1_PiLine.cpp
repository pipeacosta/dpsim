/* Copyright 2017-2021 Institute for Automation of Complex Power Systems,
 *                     EONERC, RWTH Aachen University
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 *********************************************************************************/

#include <stdexcept>

#include <dpsim-models/MathUtils.h>
#include <dpsim-models/SP/SP_Ph1_PiLine.h>

using namespace CPS;

SP::Ph1::PiLine::PiLine(String uid, String name, Logger::Level logLevel)
    : Base::Ph1::PiLine(mAttributes),
      CompositePowerComp<Complex>(uid, name, false, true, logLevel),
      mCurrent(mAttributes->create<MatrixComp>("current_vector")),
      mActivePowerBranch(mAttributes->create<Matrix>("p_branch_vector")),
      mReactivePowerBranch(mAttributes->create<Matrix>("q_branch_vector")),
      mActivePowerInjection(mAttributes->create<Real>("p_inj")),
      mReactivePowerInjection(mAttributes->create<Real>("q_inj")) {

  SPDLOG_LOGGER_INFO(mSLog, "Create {} {}", this->type(), name);
  mSLog->flush();

  // Do not create a virtual node unconditionally.
  //
  // createSubComponents() decides whether one is actually required:
  //   R > 0 && L > 0 -> one virtual node
  //   otherwise       -> direct single series element, no virtual node
  setVirtualNodeNumber(0);

  setTerminalNumber(2);

  **mIntfVoltage = MatrixComp::Zero(1, 1);
  **mIntfCurrent = MatrixComp::Zero(1, 1);
  **mCurrent = MatrixComp::Zero(2, 1);
  **mActivePowerBranch = Matrix::Zero(2, 1);
  **mReactivePowerBranch = Matrix::Zero(2, 1);
}

void SP::Ph1::PiLine::setParameters(Real resistance, Real inductance,
                                    Real capacitance, Real conductance) {

  if (resistance < 0.) {
    SPDLOG_LOGGER_ERROR(mSLog,
                        "Line {}: negative series resistance {} [Ohm] is not "
                        "supported",
                        this->name(), resistance);
    throw std::invalid_argument("PiLine series resistance must be >= 0");
  }

  if (inductance < 0.) {
    SPDLOG_LOGGER_ERROR(mSLog,
                        "Line {}: negative series inductance {} [H] is not "
                        "supported",
                        this->name(), inductance);
    throw std::invalid_argument("PiLine series inductance must be >= 0");
  }

  if (resistance == 0. && inductance == 0.) {
    SPDLOG_LOGGER_ERROR(
        mSLog,
        "Line {}: both series resistance and inductance are zero. "
        "An ideal short requires node merging / an explicit ideal-connection "
        "model and cannot be represented by PiLine MNA admittances.",
        this->name());
    throw std::invalid_argument(
        "PiLine cannot represent an ideal short with R=0 and L=0");
  }

  **mSeriesRes = resistance;
  **mSeriesInd = inductance;

  SPDLOG_LOGGER_INFO(mSLog, "Resistance={} [Ohm] Inductance={} [H]",
                     **mSeriesRes, **mSeriesInd);

  if (capacitance > 0.) {
    **mParallelCap = capacitance;
    mParallelCapIsFallback = false;
  } else {
    // Keep the existing small MNA-conditioning capacitance. The fallback is
    // excluded from standalone PF through mParallelCapIsFallback.
    **mParallelCap = 1e-12;
    mParallelCapIsFallback = true;

    SPDLOG_LOGGER_WARN(
        mSLog,
        "Zero/non-positive value for Capacitance, setting MNA fallback "
        "C={} [F]",
        **mParallelCap);
  }

  if (conductance > 0.) {
    **mParallelCond = conductance;
    mParallelCondIsFallback = false;
  } else {
    // MNA needs a small shunt conductance for conditioning. Standalone PF
    // excludes this fallback through mParallelCondIsFallback.
    **mParallelCond = 1e-6;
    mParallelCondIsFallback = true;

    SPDLOG_LOGGER_WARN(
        mSLog,
        "Zero/non-positive value for Conductance, setting MNA fallback "
        "G={} [S]",
        **mParallelCond);
  }

  SPDLOG_LOGGER_INFO(mSLog, "Capacitance={} [F] Conductance={} [S]",
                     **mParallelCap, **mParallelCond);

  mSLog->flush();
  mParametersSet = true;
}

/// DEPRECATED: Delete method
SimPowerComp<Complex>::Ptr SP::Ph1::PiLine::clone(String name) {
  auto copy = PiLine::make(name, mLogLevel);
  copy->setParameters(**mSeriesRes, **mSeriesInd, **mParallelCap,
                      **mParallelCond);
  return copy;
}

// #### Powerflow section ####

Real SP::Ph1::PiLine::getBaseVoltage() const { return mBaseVoltage; }

void SP::Ph1::PiLine::setBaseVoltage(Real baseVoltage) {
  mBaseVoltage = baseVoltage;
}

void SP::Ph1::PiLine::calculatePerUnitParameters(Real baseApparentPower,
                                                 Real baseOmega) {
  SPDLOG_LOGGER_INFO(mSLog, "#### Calculate Per Unit Parameters for {}",
                     **mName);

  mBaseApparentPower = baseApparentPower;
  mBaseOmega = baseOmega;

  SPDLOG_LOGGER_INFO(mSLog, "Base Power={} [VA]  Base Omega={} [1/s]",
                     baseApparentPower, baseOmega);

  mBaseImpedance = (mBaseVoltage * mBaseVoltage) / mBaseApparentPower;
  mBaseAdmittance = 1.0 / mBaseImpedance;
  mBaseInductance = mBaseImpedance / mBaseOmega;
  mBaseCapacitance = 1.0 / mBaseOmega / mBaseImpedance;
  mBaseCurrent =
      baseApparentPower /
      (mBaseVoltage *
       sqrt(3)); // I_base=(S_threephase/3)/(V_line_to_line/sqrt(3))

  SPDLOG_LOGGER_INFO(mSLog, "Base Voltage={} [V]  Base Impedance={} [Ohm]",
                     mBaseVoltage, mBaseImpedance);

  mSeriesResPerUnit = **mSeriesRes / mBaseImpedance;
  mSeriesIndPerUnit = **mSeriesInd / mBaseInductance;
  mParallelCapPerUnit = **mParallelCap / mBaseCapacitance;
  mParallelCondPerUnit = **mParallelCond / mBaseAdmittance;

  SPDLOG_LOGGER_INFO(mSLog, "Resistance={} [pu] Reactance={} [pu]",
                     mSeriesResPerUnit, 1. * mSeriesIndPerUnit);

  SPDLOG_LOGGER_INFO(mSLog, "Susceptance={} [pu] Conductance={} [pu]",
                     1. * mParallelCapPerUnit, mParallelCondPerUnit);

  mSLog->flush();
}

void SP::Ph1::PiLine::pfApplyAdmittanceMatrixStamp(SparseMatrixCompRow &Y) {
  int bus1 = this->matrixNodeIndex(0);
  int bus2 = this->matrixNodeIndex(1);

  // Exclude MNA-conditioning fallback shunts from standalone PF.
  Real parallelCondPerUnit =
      (mBehaviour == Behaviour::PFSimulation && mParallelCondIsFallback)
          ? 0.
          : mParallelCondPerUnit;

  Real parallelCapPerUnit =
      (mBehaviour == Behaviour::PFSimulation && mParallelCapIsFallback)
          ? 0.
          : mParallelCapPerUnit;

  // Series impedance is allowed to be purely resistive or purely inductive,
  // but setParameters() rejects the zero-zero case.
  Complex y =
      Complex(1., 0.) / Complex(mSeriesResPerUnit, 1. * mSeriesIndPerUnit);

  Complex ys =
      Complex(parallelCondPerUnit, 1. * parallelCapPerUnit) / Complex(2., 0.);

  mY_element = MatrixComp(2, 2);
  mY_element(0, 0) = y + ys;
  mY_element(0, 1) = -y;
  mY_element(1, 0) = -y;
  mY_element(1, 1) = y + ys;

  for (int i = 0; i < 2; i++) {
    for (int j = 0; j < 2; j++) {
      if (!Math::isFinite(mY_element.coeff(i, j))) {
        SPDLOG_LOGGER_ERROR(
            mSLog,
            "Line {}: non-finite per-unit admittance {} in element Y({},{})",
            this->name(), Logger::complexToString(mY_element.coeff(i, j)), i,
            j);
        throw InvalidArgumentException();
      }
    }
  }

  Y.coeffRef(bus1, bus1) += mY_element.coeff(0, 0);
  Y.coeffRef(bus1, bus2) += mY_element.coeff(0, 1);
  Y.coeffRef(bus2, bus2) += mY_element.coeff(1, 1);
  Y.coeffRef(bus2, bus1) += mY_element.coeff(1, 0);

  SPDLOG_LOGGER_INFO(mSLog, "#### PF Y matrix stamping #### ");
  SPDLOG_LOGGER_INFO(mSLog, "{}", mY_element);
  mSLog->flush();
}

void SP::Ph1::PiLine::updateBranchFlow(VectorComp &current,
                                       VectorComp &powerflow) {
  **mCurrent = current * mBaseCurrent;
  **mActivePowerBranch = powerflow.real() * mBaseApparentPower;
  **mReactivePowerBranch = powerflow.imag() * mBaseApparentPower;
}

void SP::Ph1::PiLine::storeNodalInjection(Complex powerInjection) {
  **mActivePowerInjection = std::real(powerInjection) * mBaseApparentPower;
  **mReactivePowerInjection = std::imag(powerInjection) * mBaseApparentPower;
}

MatrixComp SP::Ph1::PiLine::Y_element() { return mY_element; }

// #### MNA section ####

void SP::Ph1::PiLine::createSubComponents() {
  if (mSubCompCreated)
    return;

  mSubCompCreated = true;

  const bool hasSeriesResistance = **mSeriesRes > 0.;
  const bool hasSeriesInductance = **mSeriesInd > 0.;

  if (!hasSeriesResistance && !hasSeriesInductance) {
    throw std::invalid_argument(
        "PiLine cannot create MNA series branch with R=0 and L=0");
  }

  // Always use a finite positive MNA-conditioning shunt. This also protects
  // against a zero value reaching createSubComponents() from another setup
  // path.
  if (**mParallelCond <= 0.) {
    **mParallelCond = 1e-6;
    mParallelCondIsFallback = true;
  }

  // -----------------------------------------------------------------------
  // Series branch
  // -----------------------------------------------------------------------
  //
  // R > 0, L > 0:
  //
  //   terminal 0 --- R --- vnode --- L --- terminal 1
  //
  // Exactly one of R/L is zero:
  //
  //   terminal 0 --- remaining element --- terminal 1
  //
  // This is the key fix: never instantiate a zero-ohm resistor or a
  // zero-henry inductor merely to preserve the composite topology.
  if (hasSeriesResistance && hasSeriesInductance) {
    setVirtualNodeNumber(1);

    mSubSeriesResistor =
        std::make_shared<SP::Ph1::Resistor>(**mName + "_res", mLogLevel);
    mSubSeriesResistor->setParameters(**mSeriesRes);
    mSubSeriesResistor->connect({mTerminals[0]->node(), mVirtualNodes[0]});

    addMNASubComponent(mSubSeriesResistor, MNA_SUBCOMP_TASK_ORDER::NO_TASK,
                       MNA_SUBCOMP_TASK_ORDER::TASK_BEFORE_PARENT, false);

    mSubSeriesInductor =
        std::make_shared<SP::Ph1::Inductor>(**mName + "_ind", mLogLevel);
    mSubSeriesInductor->setParameters(**mSeriesInd);
    mSubSeriesInductor->connect({mVirtualNodes[0], mTerminals[1]->node()});

    addMNASubComponent(mSubSeriesInductor, MNA_SUBCOMP_TASK_ORDER::NO_TASK,
                       MNA_SUBCOMP_TASK_ORDER::TASK_BEFORE_PARENT, true);

  } else if (hasSeriesInductance) {
    // Purely inductive series branch (R = 0).
    // No virtual node and, importantly, no zero-ohm resistor.
    mSubSeriesInductor =
        std::make_shared<SP::Ph1::Inductor>(**mName + "_ind", mLogLevel);
    mSubSeriesInductor->setParameters(**mSeriesInd);
    mSubSeriesInductor->connect(
        {mTerminals[0]->node(), mTerminals[1]->node()});

    addMNASubComponent(mSubSeriesInductor, MNA_SUBCOMP_TASK_ORDER::NO_TASK,
                       MNA_SUBCOMP_TASK_ORDER::TASK_BEFORE_PARENT, true);

  } else {
    // Purely resistive series branch (L = 0).
    // No virtual node and no zero-henry inductor.
    mSubSeriesResistor =
        std::make_shared<SP::Ph1::Resistor>(**mName + "_res", mLogLevel);
    mSubSeriesResistor->setParameters(**mSeriesRes);
    mSubSeriesResistor->connect(
        {mTerminals[0]->node(), mTerminals[1]->node()});

    addMNASubComponent(mSubSeriesResistor, MNA_SUBCOMP_TASK_ORDER::NO_TASK,
                       MNA_SUBCOMP_TASK_ORDER::TASK_BEFORE_PARENT, false);
  }

  // -----------------------------------------------------------------------
  // Parallel shunts
  // -----------------------------------------------------------------------

  mSubParallelResistor0 =
      std::make_shared<SP::Ph1::Resistor>(**mName + "_con0", mLogLevel);
  mSubParallelResistor0->setParameters(2. / **mParallelCond);
  mSubParallelResistor0->connect(
      SimNode::List{SimNode::GND, mTerminals[0]->node()});

  addMNASubComponent(mSubParallelResistor0, MNA_SUBCOMP_TASK_ORDER::NO_TASK,
                     MNA_SUBCOMP_TASK_ORDER::TASK_BEFORE_PARENT, false);

  mSubParallelResistor1 =
      std::make_shared<SP::Ph1::Resistor>(**mName + "_con1", mLogLevel);
  mSubParallelResistor1->setParameters(2. / **mParallelCond);
  mSubParallelResistor1->connect(
      SimNode::List{SimNode::GND, mTerminals[1]->node()});

  addMNASubComponent(mSubParallelResistor1, MNA_SUBCOMP_TASK_ORDER::NO_TASK,
                     MNA_SUBCOMP_TASK_ORDER::TASK_BEFORE_PARENT, false);

  if (**mParallelCap > 0.) {
    mSubParallelCapacitor0 =
        std::make_shared<SP::Ph1::Capacitor>(**mName + "_cap0", mLogLevel);
    mSubParallelCapacitor0->setParameters(**mParallelCap / 2.);
    mSubParallelCapacitor0->connect(
        SimNode::List{SimNode::GND, mTerminals[0]->node()});

    addMNASubComponent(mSubParallelCapacitor0,
                       MNA_SUBCOMP_TASK_ORDER::NO_TASK,
                       MNA_SUBCOMP_TASK_ORDER::TASK_BEFORE_PARENT, true);

    mSubParallelCapacitor1 =
        std::make_shared<SP::Ph1::Capacitor>(**mName + "_cap1", mLogLevel);
    mSubParallelCapacitor1->setParameters(**mParallelCap / 2.);
    mSubParallelCapacitor1->connect(
        SimNode::List{SimNode::GND, mTerminals[1]->node()});

    addMNASubComponent(mSubParallelCapacitor1,
                       MNA_SUBCOMP_TASK_ORDER::NO_TASK,
                       MNA_SUBCOMP_TASK_ORDER::TASK_BEFORE_PARENT, true);
  }
}

void SP::Ph1::PiLine::initializeParentFromNodesAndTerminals(Real frequency) {
  const Real omega = 2. * PI * frequency;
  const Complex impedance = {**mSeriesRes, omega * **mSeriesInd};

  if (std::abs(impedance) == 0.) {
    throw std::runtime_error(
        "PiLine initialization encountered zero series impedance");
  }

  (**mIntfVoltage)(0, 0) =
      initialSingleVoltage(1) - initialSingleVoltage(0);

  (**mIntfCurrent)(0, 0) = (**mIntfVoltage)(0, 0) / impedance;

  // A virtual node exists only for the R+L representation.
  if (hasVirtualNodes()) {
    mVirtualNodes[0]->setInitialVoltage(initialSingleVoltage(0) +
                                        (**mIntfCurrent)(0, 0) * **mSeriesRes);

    SPDLOG_LOGGER_INFO(
        mSLog,
        "\n--- Initialization from powerflow ---"
        "\nVoltage across: {:s}"
        "\nCurrent: {:s}"
        "\nTerminal 0 voltage: {:s}"
        "\nTerminal 1 voltage: {:s}"
        "\nVirtual Node 1 voltage: {:s}"
        "\n--- Initialization from powerflow finished ---",
        Logger::phasorToString((**mIntfVoltage)(0, 0)),
        Logger::phasorToString((**mIntfCurrent)(0, 0)),
        Logger::phasorToString(initialSingleVoltage(0)),
        Logger::phasorToString(initialSingleVoltage(1)),
        Logger::phasorToString(mVirtualNodes[0]->initialSingleVoltage()));
  } else {
    SPDLOG_LOGGER_INFO(
        mSLog,
        "\n--- Initialization from powerflow ---"
        "\nVoltage across: {:s}"
        "\nCurrent: {:s}"
        "\nTerminal 0 voltage: {:s}"
        "\nTerminal 1 voltage: {:s}"
        "\nSeries representation: direct {}"
        "\n--- Initialization from powerflow finished ---",
        Logger::phasorToString((**mIntfVoltage)(0, 0)),
        Logger::phasorToString((**mIntfCurrent)(0, 0)),
        Logger::phasorToString(initialSingleVoltage(0)),
        Logger::phasorToString(initialSingleVoltage(1)),
        mSubSeriesInductor ? "inductor" : "resistor");
  }
}

void SP::Ph1::PiLine::mnaParentAddPostStepDependencies(
    AttributeBase::List &prevStepDependencies,
    AttributeBase::List &attributeDependencies,
    AttributeBase::List &modifiedAttributes,
    Attribute<Matrix>::Ptr &leftVector) {

  attributeDependencies.push_back(leftVector);
  modifiedAttributes.push_back(mIntfVoltage);
  modifiedAttributes.push_back(mIntfCurrent);
}

void SP::Ph1::PiLine::mnaParentPostStep(
    Real time, Int timeStepCount, Attribute<Matrix>::Ptr &leftVector) {

  this->mnaUpdateVoltage(**leftVector);
  this->mnaUpdateCurrent(**leftVector);
}

void SP::Ph1::PiLine::mnaCompUpdateVoltage(const Matrix &leftVector) {
  (**mIntfVoltage)(0, 0) = 0;

  if (terminalNotGrounded(1))
    (**mIntfVoltage)(0, 0) =
        Math::complexFromVectorElement(leftVector, matrixNodeIndex(1));

  if (terminalNotGrounded(0))
    (**mIntfVoltage)(0, 0) =
        (**mIntfVoltage)(0, 0) -
        Math::complexFromVectorElement(leftVector, matrixNodeIndex(0));
}

void SP::Ph1::PiLine::mnaCompUpdateCurrent(const Matrix &leftVector) {
  if (mSubSeriesInductor) {
    (**mIntfCurrent)(0, 0) = mSubSeriesInductor->intfCurrent()(0, 0);
    return;
  }

  if (mSubSeriesResistor) {
    (**mIntfCurrent)(0, 0) = mSubSeriesResistor->intfCurrent()(0, 0);
    return;
  }

  throw std::runtime_error("PiLine has no valid series subcomponent");
}

// #### Tear methods ####

MNAInterface::List SP::Ph1::PiLine::mnaTearGroundComponents() {
  MNAInterface::List gndComponents;

  if (mSubParallelResistor0)
    gndComponents.push_back(mSubParallelResistor0);

  if (mSubParallelResistor1)
    gndComponents.push_back(mSubParallelResistor1);

  if (mSubParallelCapacitor0)
    gndComponents.push_back(mSubParallelCapacitor0);

  if (mSubParallelCapacitor1)
    gndComponents.push_back(mSubParallelCapacitor1);

  return gndComponents;
}

void SP::Ph1::PiLine::mnaTearInitialize(Real omega, Real timeStep) {
  if (mSubSeriesResistor) {
    mSubSeriesResistor->mnaTearSetIdx(mTearIdx);
    mSubSeriesResistor->mnaTearInitialize(omega, timeStep);
  }

  if (mSubSeriesInductor) {
    mSubSeriesInductor->mnaTearSetIdx(mTearIdx);
    mSubSeriesInductor->mnaTearInitialize(omega, timeStep);
  }
}

void SP::Ph1::PiLine::mnaTearApplyMatrixStamp(SparseMatrixRow &tearMatrix) {
  if (mSubSeriesResistor)
    mSubSeriesResistor->mnaTearApplyMatrixStamp(tearMatrix);

  if (mSubSeriesInductor)
    mSubSeriesInductor->mnaTearApplyMatrixStamp(tearMatrix);
}

void SP::Ph1::PiLine::mnaTearApplyVoltageStamp(Matrix &voltageVector) {
  if (mSubSeriesInductor)
    mSubSeriesInductor->mnaTearApplyVoltageStamp(voltageVector);
}

void SP::Ph1::PiLine::mnaTearPostStep(Complex voltage, Complex current) {
  if (mSubSeriesInductor) {
    mSubSeriesInductor->mnaTearPostStep(voltage - current * **mSeriesRes,
                                        current);
  }

  if (mSubSeriesResistor) {
    mSubSeriesResistor->mnaTearPostStep(voltage, current);
  }

  // Keep the parent interface variables consistent even though the SP
  // resistor/inductor tear post-step implementations are no-ops.
  (**mIntfVoltage)(0, 0) = voltage;
  (**mIntfCurrent)(0, 0) = current;
}
