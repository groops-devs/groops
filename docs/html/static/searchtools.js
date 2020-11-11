/** Global variables */

var idx = lunr(function () {
              this.ref('key')
              this.field('name')
              this.field('display_text')
              this.field('description')
              this.field('config_table')
              this.metadataWhitelist = ['position']

              for(var i in documents)
                this.add(documents[i]);
            })

var wrapThreshold = 6  /** number of search results per page */
var paginationCount = 5  /** number of (maximum) pagination links */

var activePageIndex = 0  /** index of currently viewed search result page */
var currentPaginationLinks;  /** array of currently available page links in pagination list */
var resultRangePerPage;  /** array containing the intervals of search results per page (ctd. blockIndex) */

/** Functions for search result presentation */

/**
 * Show/hide search result boxes
 *
 * @param oldIndex page index of results to hide
 * @param newIndex page index of results to show
 */
function toggleSearchResults(oldIndex, newIndex) {

    var resultList = document.getElementById("searchResults");
    for(var k = resultRangePerPage[oldIndex]; k<resultRangePerPage[oldIndex+1]; k++)
        resultList.childNodes[k].style.display = "none";

    for(var k = resultRangePerPage[newIndex]; k<resultRangePerPage[newIndex+1]; k++)
        resultList.childNodes[k].style.display = "";
}

/**
 * Navigation of search results with pagination buttons.
 * Shows/hides search result boxes and updates links and meta-data
 *
 * @param resultRangePerPage  array containing the intervals of search results per page (ctd. blockIndex)
 * @param pageIndex index of page to be displayed
 */
function searchResultNavigation(resultRangePerPage, pageIndex) {

    var newPageIndex = Math.min(resultRangePerPage.length-2, Math.max(0, pageIndex))  // make sure we are in the correct range
    if(newPageIndex == activePageIndex)  // short curcuit if the same page is requested
        return;

    toggleSearchResults(activePageIndex, newPageIndex)  // hide/show the requested search results

    var currentCenter = currentPaginationLinks[Math.floor(currentPaginationLinks.length/2)]  // center of current pagination bar
    var shift = newPageIndex - currentCenter  // new center is requested page

    if( (currentPaginationLinks[currentPaginationLinks.length-1]+shift)>(resultRangePerPage.length-2) || (currentPaginationLinks[0]+shift)<0 )  // if the shift would put outside bounds, ignore
        shift = 0

    for(i in currentPaginationLinks)  // update pagination link numbers
        currentPaginationLinks[i] += shift

    activePageIndex = newPageIndex  // update HTML elements
    updatePaginationLinks(currentPaginationLinks, activePageIndex)
}

/**
 * Update HTML elements of pagination buttons
 *
 * @param paginationLinks  array of page indices currently in pagination list
 * @param newIndex index of page to be displayed
 */
function updatePaginationLinks(paginationLinks, newIndex) {

    var paginationSpan = document.getElementById("pagination-page-number")
    paginationSpan.innerHTML = (newIndex+1)

    for(var k = 0; k<paginationLinks.length; k++)
    {
        var item = document.getElementById("pagination-link-"+k)
        item.setAttribute("onclick", "searchResultNavigation(resultRangePerPage, "+paginationLinks[k]+");")
        item.innerHTML = '<span class="'+(paginationLinks[k] == newIndex ? 'text-warning' :'text-primary')+' page-link">'+(paginationLinks[k]+1)+'</a>'
    }

}

/**
 * Emphasize text based on search result position
 *
 * @param position array of search result position in text
 * @param text text to be changed
 */
function highlightResult(position, text) {

    var highlight = '<em class="text-warning">' + text.slice(position[0], position[0]+position[1]+1) + '</em>'
    return text.slice(0, position[0]) + highlight + text.slice(position[0]+position[1]+1, text.length)
}

/**
 * Function to expand a search result box on button click
 *
 * @param hiddedId id of box to be expanded
 */
function expandSearchBox(hiddenId)
{
    var searchBox = document.getElementById(hiddenId)
    searchBox.classList.toggle("faded")
}

/**
 * Parse search terms for URL. This function generates the query which is passed to the search module
 *
 * @param rawInput URL POST token
 */
function parseSearchTerms(rawInput)
{
    var tokens = rawInput.split('+')

    var searchString = tokens[0]
    for(var i = 1; i<tokens.length; i++)
        searchString += ' '+tokens[i]

    for(var i = 0; i<tokens.length; i++)
        searchString += '  name:*'+tokens[i]+'*'

    return searchString
}

/** Mustache templates for reoccurring HTML elements */

/** HTML template for a search result box */
var templateSearchResult = '<div class="card border-light mb-3" style="{{style}}"><div class="card-header">' +
                           '<a href="{{key}}.html">{{name}}</a>' +
                           '<span class="badge badge-primary float-right" onclick="expandSearchBox(\'search:{{name}}\');">+</span></div>' +
                           '<div id="search:{{name}}" class="card-body search-result faded"><p>{{{content}}}</p></div></div>'

/** HTML template for the pagination links */
var templatePagination = '<ul class="pagination">'+
                         '<li class="page-item" onclick="searchResultNavigation(resultRangePerPage, activePageIndex-1);"><span class="text-primary page-link">&laquo;</span></li>'+
                         '{{#pages}}<li class="page-item" id="pagination-link-{{.}}">{{.}}</li>{{/pages}}'+
                         '<li class="page-item" onclick="searchResultNavigation(resultRangePerPage, activePageIndex+1);"><span class="text-primary page-link">&raquo;</span></li>'+
                         '</ul>'

/** Search result generation */

/**
 * This function should be called when the search page is loaded. It parses the POST data from the URL and queries the
 * search index. The results are dynamically displayed as (potentially) hidden cards.
 */
function startSearch() {

    var rawSearchString = window.location.search.substr(1).split('=')[1]

    var metaInfo = document.getElementById("search-meta-info")
    var list = document.getElementById("searchResults")
    while(list.firstChild) {
        list.removeChild(list.firstChild);
    }

    if(rawSearchString == '')  // short circuit on empty search
    {
        metaInfo.innerHTML = 'Empty search.';
        return;
    }

    var searchQuery = parseSearchTerms(rawSearchString)  // create query from user input and search index
    var results = idx.search(searchQuery);

    if(results.length == 0)  // short circuit on no results
    {
        metaInfo.innerHTML = 'Your search - <b>"'+rawSearchString+'"</b> - did not return any results.'
        return;
    }

    resultRangePerPage = [0]
    while(resultRangePerPage[resultRangePerPage.length-1]<results.length)
        resultRangePerPage.push(Math.min(resultRangePerPage[resultRangePerPage.length-1]+wrapThreshold, results.length))

    metaInfo.innerHTML = 'Page <span id="pagination-page-number">1</span> of '+(resultRangePerPage.length-1)+' ('+results.length+' results).'

    for(i in results)
    {
        var metadata = results[i].matchData.metadata;
        var positions =  {'name': [], 'description': [], 'display_text': [], 'config_table': []}
        for(j in metadata)
            for(k in metadata[j])
               positions[k] = positions[k].concat(metadata[j][k].position);

        var fullText = documents[results[i].ref].display_text

        var data = {'name': documents[results[i].ref].name, 'key': results[i].ref,'content':  fullText/*documents[results[i].ref].description*/, 'style': (i>=wrapThreshold) ? 'display: none;' : ''}
        var renderedSearchResult = Mustache.render(templateSearchResult, data)

        list.insertAdjacentHTML('beforeend', renderedSearchResult)
    }

    currentPaginationLinks = Array.apply(null, {length: Math.min(paginationCount, resultRangePerPage.length-1)}).map(Number.call, Number)
    var data = {'pages': currentPaginationLinks}
    var renderedPagination = Mustache.render(templatePagination, data)

    document.getElementById("content").insertAdjacentHTML('beforeend', renderedPagination)

    activePageIndex = 0
    updatePaginationLinks(currentPaginationLinks, activePageIndex)
 }
