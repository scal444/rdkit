//
// Copyright (C) David Cosgrove 2025.
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#ifndef CONTROLCHANDLER_H
#define CONTROLCHANDLER_H

#include <atomic>
#include <csignal>
#include <mutex>
#include <stdexcept>

namespace RDKit {

class ControlCCaught : public std::runtime_error {
 public:
  explicit ControlCCaught()
      : std::runtime_error("The process was interrupted with Ctrl+c") {};
};

//! This class catches a control-C/SIGINT and sets the flag d_gotSignal
//! if one is received.  It is intended to be used inside a long
//! C++ calculation called from Python which intercepts the signal
//! handler.  The C++ code must check the value of d_gotSignal
//! periodically and act accordingly.  The destructor resets
//! the signal handler and flag after the last concurrent user exits.
//! Example usage, inside a boost::python wrapper:
//!  ResultsObject results;
//!  {
//!   NOGIL gil;
//!    results = someFunction();
//!  }
//!  if (results.getCancelled()) {
//!    throw_runtime_error("someFunction cancelled");
//!  }
//! It's important that the exception is thrown once the GIL has been
//! released, otherwise a crash is inevitable at some future point.
class ControlCHandler {
 public:
  ControlCHandler() {
    std::lock_guard<std::mutex> lock(d_mutex);
    if (d_handlerCount++ == 0) {
      d_gotSignal.clear(std::memory_order_relaxed);
      d_prev_handler = std::signal(SIGINT, signalHandler);
    }
  }
  ControlCHandler(const ControlCHandler &) = delete;
  ControlCHandler(ControlCHandler &&) = delete;
  ControlCHandler &operator=(const ControlCHandler &) = delete;
  ControlCHandler &operator=(ControlCHandler &&) = delete;
  ~ControlCHandler() {
    std::lock_guard<std::mutex> lock(d_mutex);
    if (--d_handlerCount == 0) {
      std::signal(SIGINT, d_prev_handler);
      d_gotSignal.clear(std::memory_order_relaxed);
    }
  }
  static bool getGotSignal() {
    return d_gotSignal.test(std::memory_order_relaxed);
  }
  static void signalHandler(int signalNumber) {
    if (signalNumber == SIGINT) {
      d_gotSignal.test_and_set(std::memory_order_relaxed);
    }
  }
  static void reset() { d_gotSignal.clear(std::memory_order_relaxed); }

 private:
  inline static std::atomic_flag d_gotSignal{};
  inline static std::mutex d_mutex;
  inline static unsigned int d_handlerCount{0};
  inline static void (*d_prev_handler)(int){SIG_DFL};
};
}  // namespace RDKit
#endif  // CONTROLCHANDLER_H
