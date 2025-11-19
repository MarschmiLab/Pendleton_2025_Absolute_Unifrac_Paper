library(ggplot2)
library(NatParksPalettes)

ggpreview <- function(...) {
  fname <- tempfile(fileext = ".png")
  ggsave(filename = fname, ...)
  system2("open", fname)
  invisible(NULL)
}

ont_theme <- theme_classic() + 
  theme(axis.text = element_text(size = 10),
        axis.title = element_text(size=12),
        strip.text = element_text(size = 10),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 12),
        strip.background = element_rect(colour = NA, fill = 'transparent'),
        panel.background = element_rect(fill='transparent'), #transparent panel bg
        plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
        panel.grid.major = element_blank(), #remove major gridlines
        panel.grid.minor = element_blank(), #remove minor gridlines
        legend.background = element_rect(fill='transparent'), #transparent legend bg
        legend.box.background = element_rect(fill='transparent') #transparent legend panel
  )

theme_set(ont_theme)

triglav <- natparks.pals(name = "Triglav", n = 6)

volcano <- rev(natparks.pals(name = "Volcanoes", n = 7))[-1]
volcano
