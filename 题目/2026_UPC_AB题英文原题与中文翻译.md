# 2026 UPC AB题英文原题与中文翻译（只保留B题）

每题先列英文原题，随后给出中文译文。中文用于辅助理解，题意以英文原文为准。

## Problem B: Predicting FIFA Penalty Trajectories

### 英文原题

In the game of FIFA football/soccer, when a penalty kick is taken, one player takes a single shot at the goal, which is only defended by the opposing team's goalkeeper. The kick is taken with the ball at the penalty mark, a point centered 11 meters from the goal, which is a rectangle 7.32 meters wide and 2.44 meters tall. Working in meters, let us take the penalty mark as our origin (0,0,0) and let the goal be in the x = 11 plane with -3.66 < y < 3.66 and 0 < z < 2.44.

To block a kick the goalkeeper must quickly anticipate the trajectory of the ball. Let the kick take place at t = 0 s when the center of the ball is at (0, 0, 0.11). Develop and test a general method for predicting the point at which the ball is most likely to enter the goal given the coordinates of the ball at t = 0.05 s and t=0.10 s, and estimating the uncertainty of your prediction.

Apply your method to a scenario where at time t = 0.05 s the ball is at (1.50, 0.45, 0.44) and at time t = 0.10 s the ball is at (2.96, 0.86, 0.74).

### 中文翻译

**B题 FIFA点球轨迹预测**

在国际足联（FIFA）足球比赛中，罚点球时，一名球员向球门进行一次射门，而球门仅由对方球队的守门员防守。射门时，足球放置在罚球点（penalty mark）；该点位于球门正中前方 11 m 处。球门为宽 7.32 m、高 2.44 m 的矩形。以米为长度单位，将罚球点取为坐标原点 (0,0,0)，并令球门位于 x = 11 的平面内，其范围为 -3.66 < y < 3.66 且 0 < z < 2.44。

为扑救点球，守门员必须迅速预判球的轨迹。设踢球发生在 t = 0 s，此时球心位于 (0, 0, 0.11)。建立并检验一种通用方法：根据球在 t = 0.05 s 和 t = 0.10 s 时的位置坐标，预测球最有可能进入球门的位置，并估计该预测的不确定性（uncertainty）。

将你的方法应用于以下情形：在 t = 0.05 s 时，球位于 (1.50, 0.45, 0.44)；在 t = 0.10 s 时，球位于 (2.96, 0.86, 0.74)。

## 来源

[University Physics Competition 2026 Competition](https://university-physics-competition.com/2026-competition)

中文译文为 AI 辅助翻译。整理日期：2026年11月7日。
