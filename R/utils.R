setTMPEnvVt <- function(inYamlFile){
  inEnvs <- yaml::read_yaml(inYamlFile)
  withr::with_envvar(new=c("vt" =inEnvs$vt), Sys.getenv("vt"))
}
