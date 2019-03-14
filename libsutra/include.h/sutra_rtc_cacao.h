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

template <typename Tin, typename Tcomp, typename Tout>
class sutra_rtc_cacao : public sutra_rtc<Tin, Tcomp, Tout> {
 private:
  std::shared_ptr<ipc::Cacao<Tin>> iRawFrame_;
  std::shared_ptr<ipc::Cacao<float>> iCalFrame_;
  std::shared_ptr<ipc::Cacao<Tcomp>> iLoopFrame_;
  std::shared_ptr<ipc::Cacao<Tout>> iCommands_;

  std::string iRawFrame_name_;
  std::string iCalFrame_name_;
  std::string iLoopFrame_name_;
  std::string iCommands_name_;

  long framecounter_;

  int nslp_;
  int ncmd_;
  int nvalid_;

  bool is_initialised_;

 public:
  sutra_rtc_cacao(std::string iCalFrame_name, std::string iLoopFrame_name);
  ~sutra_rtc_cacao();

  void publish();

 private:
  void allocateBuffers();
};

#endif /* SUTRA_RTC_CACAO_H_ */
