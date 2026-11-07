/* Independently implemented ideal-mixture atmosphere helper, SI throughout. */
function airProperties(temperatureC, pressurePa, relativeHumidityPercent) {
  'use strict';
  const values=[temperatureC,pressurePa,relativeHumidityPercent];
  if(!values.every(Number.isFinite)||temperatureC<0||temperatureC>40||pressurePa<50000||pressurePa>110000||relativeHumidityPercent<0||relativeHumidityPercent>100)
    throw Error('物性助手范围：0–40 °C，外界绝对压力 50000–110000 Pa，相对湿度 0–100%。超出范围请直接提供有来源的ρ、μ。');
  const T=temperatureC+273.15;
  const saturationPressurePa=611.2*Math.exp(17.67*temperatureC/(243.5+temperatureC));
  const vaporPressurePa=relativeHumidityPercent*.01*saturationPressurePa;
  const rho=(pressurePa-vaporPressurePa)/(287.05*T)+vaporPressurePa/(461.5*T);
  const mu=1.716e-5*Math.pow(T/273,1.5)*(273+111)/(T+111);
  return {rho,mu,saturationPressurePa,vaporPressurePa,temperatureC,pressurePa,relativeHumidityPercent,
    method:'液水平衡饱和蒸汽压近似 + 干空气/水蒸气理想混合密度；黏度采用干空气Sutherland近似，未加入湿度修正',
    source:['https://essd.copernicus.org/articles/13/155/2021/','https://doc.comsol.com/6.3/doc/com.comsol.help.cfd/cfd_ug_fluidflow_high_mach.08.43.html'],
    assumptions:['现场气象未知，输入为情景','pressurePa是环境绝对压力，不是球内充气压力','不是CIPM计量级湿空气模型','点击应用后覆盖直接输入ρ、μ；湿度不表示球皮湿润']};
}
if(typeof module!=='undefined') module.exports={airProperties};
