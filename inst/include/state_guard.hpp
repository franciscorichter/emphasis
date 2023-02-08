#ifndef EMPHASIS_GUARDS_HPP_INCLUDED
#define EMPHASIS_GUARDS_HPP_INCLUDED

#include "model.hpp"

namespace emphasis {


  class state_guard
  {
  public:
    state_guard(const state_guard&) = delete;
    state_guard& operator=(const state_guard&) = delete;

    state_guard(const Model* model)
      : model_(model), state_(nullptr)
    {}

    ~state_guard() { model_->free_state(&state_); }
    operator void** () noexcept { return &state_; }
    void invalidate_state() { model_->invalidate_state(&state_); }

  private:
    const Model* model_;
    void* state_;
  };
}

#endif
