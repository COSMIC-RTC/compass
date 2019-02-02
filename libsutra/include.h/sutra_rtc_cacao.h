/*
 * sutra_rtc_cacao.h
 *
 *  Created on: Feb 1, 2019
 *      Author: sevin
 */

#ifndef SUTRA_RTC_CACAO_H_
#define SUTRA_RTC_CACAO_H_

#include <Cacao.h>

#include <sutra_rtc.h>
#include <sutra_target.h>
#include <sutra_wfs.h>

template <typename T>
class sutra_rtc_cacao : public sutra_rtc {
 private:
  std::shared_ptr<ipc::Cacao<T>> iCalFrame_;
  std::shared_ptr<ipc::Cacao<T>> iLoopFrame_;

  std::string iCalFrame_name_;
  std::string iLoopFrame_name_;

  long framecounter_;

  int nslp_;
  int ncmd_;
  int nvalid_;

  bool is_initialised_;

 public:
  sutra_rtc_cacao(carma_context *context, std::string iCalFrame_name,
                  std::string iLoopFrame_name);
  ~sutra_rtc_cacao();

  void publish();

 private:
  void allocateBuffers();
};

#endif /* SUTRA_RTC_CACAO_H_ */
