library(tidyverse)


# You need to load a package (like magrittr or dplyr) that defines the function first, then it should work.

install.packages("magrittr") # package installations are only needed the first time you use it
install.packages("dplyr")    # alternative installation of the %>%
library(magrittr) # needs to be run every time you start R and want to use %>%
library(dplyr)    # alternatively, this also loads %>%

# gapminder data in wide form from Carpentries
gapminder_w_url <- "https://bit.ly/2vEDq5b"
# read the data into dataframe
gapminder_wide <- read_csv(gapminder_w_url)
# first3 rows
head(gapminder_wide, n=3)

## # A tibble: 3 x 38
##   continent country gdpPercap_1952 gdpPercap_1957 gdpPercap_1962
##   <chr>     <chr>            <dbl>          <dbl>          <dbl>
## 1 Africa    Algeria          2449.          3014.          2551.
## 2 Africa    Angola           3521.          3828.          4269.
## 3 Africa    Benin            1063.           960.           949.
## # … with 33 more variables: gdpPercap_1967 <dbl>, gdpPercap_1972 <dbl

gapminder_life <- gapminder_wide %>%
  filter(continent %in% c("Africa","Europe")) %>%
  select(continent,country,starts_with('lifeExp'))

gapminder_life <- gapminder_life %>% 
  unite("continent_country", c(continent,country)) %>%
  column_to_rownames("continent_country")
head(gapminder_life)